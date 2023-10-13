//
// Created by jonah on 3/26/21.
//


#include "CUDADNANMInteraction.h"
#include "CUDA_DNANM.cuh"
//#include "CUDA_DNA.cuh"
//#include "../Lists/CUDASimpleVerletList.h"
//#include "../Lists/CUDANoList.h"
#include "../../Interactions/DNANMInteraction.h"


CUDADNANMInteraction::CUDADNANMInteraction(bool btp) : DNANMInteraction(btp), CUDABaseInteraction() {
    
    _read_par = false;
    _angular = btp;
    //Not copied over to device memory
    _spring_potential = NULL;
    _spring_eqdist = NULL;
    _affected_len = NULL;

    //Copied over to device memory
    _h_ang_params = NULL;
    _d_ang_params = NULL;

    _h_ang_kbkt = NULL;
    _d_ang_kbkt = NULL;

    _h_affected_indx = NULL;
    _d_affected_indx = NULL;

    _h_affected = NULL;
    _d_affected = NULL;

    _h_aff_eqdist = NULL;
    _d_aff_eqdist = NULL;

    _h_aff_gamma = NULL;
    _d_aff_gamma = NULL;

    _spring_param_size_number = 0;
    _ang_param_size = 0;

    _edge_compatible = true;
    _use_debye_huckel = false;
    _use_oxDNA2_coaxial_stacking = false;
    _use_oxDNA2_FENE = false;
}


CUDADNANMInteraction::~CUDADNANMInteraction() {
    //Delete All pointers required for spring potential parameters
    if(_spring_potential != NULL) delete[] _spring_potential;
    if(_spring_eqdist != NULL) delete[] _spring_eqdist;

    if(_h_ang_params != NULL) delete[] _h_ang_params;
    if(_h_ang_kbkt != NULL) delete[] _h_ang_kbkt;
    if(_d_ang_params != NULL) CUDA_SAFE_CALL( cudaFree(_d_ang_params) );
    if(_d_ang_kbkt != NULL) CUDA_SAFE_CALL( cudaFree(_d_ang_kbkt) );

    if(_affected_len != NULL) delete[] _affected_len;

    if(_d_affected != NULL) CUDA_SAFE_CALL(cudaFree(_d_affected));
    if(_d_aff_gamma != NULL) CUDA_SAFE_CALL(cudaFree(_d_aff_gamma));
    if(_d_aff_eqdist != NULL) CUDA_SAFE_CALL(cudaFree(_d_aff_eqdist));
    if(_d_affected_indx != NULL) CUDA_SAFE_CALL( cudaFree(_d_affected_indx) );

    if(_h_affected != NULL) delete[] _h_affected;
    if(_h_aff_gamma != NULL) delete[] _h_aff_gamma;
    if(_h_aff_eqdist != NULL) delete[] _h_aff_eqdist;
    if(_h_affected_indx != NULL) delete[] _h_affected_indx;

    if(_d_is_strand_end != nullptr) {
        CUDA_SAFE_CALL(cudaFree(_d_is_strand_end));
    }
}


void CUDADNANMInteraction::get_settings(input_file &inp) {
    std::string inter_type;
    if (getInputString(&inp, "parfile", _parameterfile, 0) != KEY_FOUND) {
        throw oxDNAException("Key 'parfile' not found. Necessary for Protein sims.");
    }

    char s[5] = "none";
    if(strcmp(this->_parameterfile, s) != 0) _read_par = true;

    if (!getInputString(&inp, "topology", this->_topology_filename, 0) == KEY_FOUND){
        throw oxDNAException("Key 'topology_file' not found.");
    }
    DNANMInteraction::get_settings(inp);

    _use_debye_huckel = true;
    _use_oxDNA2_coaxial_stacking = true;
    _use_oxDNA2_FENE = true;

    // we don't need the F4_... terms as the macros are used in the CUDA_DNA.cuh file; this doesn't apply for the F2_K term
    F2_K[1] = CXST_K_OXDNA2;
    _debye_huckel_half_charged_ends = true;
    this->_grooving = true;
    // end copy from DNA2Interaction

    // copied from DNA2Interaction::get_settings() (CPU), the least bad way of doing things
    getInputNumber(&inp, "salt_concentration", &_salt_concentration, 1);
    getInputBool(&inp, "dh_half_charged_ends", &_debye_huckel_half_charged_ends, 0);

    // lambda-factor (the dh length at T = 300K, I = 1.0)
    _debye_huckel_lambdafactor = 0.3616455f;
    getInputFloat(&inp, "dh_lambda", &_debye_huckel_lambdafactor, 0);

    // the prefactor to the Debye-Huckel term
    _debye_huckel_prefactor = 0.0543f;
    getInputFloat(&inp, "dh_strength", &_debye_huckel_prefactor, 0);
}


void CUDADNANMInteraction::cuda_init(int N) {
    CUDABaseInteraction::cuda_init(N);
    DNANMInteraction::init();

//    Addition of Reading Parameter File -> Moved from get_settings due to needing to fill variables that are filled in the CPU version of DNANMInteraction::read_topology
    std::fstream top;
    int tmp1, tmp2;
    top.open(this->_topology_filename, std::ios::in);
    if (top.is_open()){
        top >> tmp1 >> tmp2 >> this->ndna >> this->npro >> this->ndnas;
        top >> this->_firststrand;
        top.close();
    } else {
        throw oxDNAException("Could not open Topology File");
    }


    if(this->_firststrand < 0) offset = 0;
    else if(this->_firststrand > 0) offset = this->ndna;
    else throw oxDNAException("No Strand should have an ID of 0");


    if(_read_par){
        //Initalizing Some Host and Device Arrays for Spring Parameters
        _spring_param_size_number = sizeof(c_number) * (this->npro*this->npro);
        _ang_param_size = sizeof(c_number) * (this->npro*4);

        _spring_potential = new c_number[this->npro*this->npro]();
        _spring_eqdist = new c_number[this->npro*this->npro]();

        if(_angular) { // only need for anmt
            //Initializing Host and Device Arrays for Angular Parameters
            _h_ang_params = new c_number[this->npro*4]();
            CUDA_SAFE_CALL( cudaMalloc(&_d_ang_params, _ang_param_size));
            if (this->_parameter_kbkt) {
                _h_ang_kbkt = new c_number[this->npro * 4]();
                CUDA_SAFE_CALL(cudaMalloc(&_d_ang_kbkt, sizeof(c_number) * (this->npro * 2)));
            }
        }

        char potswitch = 'x';
        c_number potential = 0.f;
        c_number dist = 0.f;
        //initializing array members
        for(int i = 0; i< (this->npro*this->npro); i++){
            _spring_eqdist[i] = dist;
            _spring_potential[i] = potential;
            if(_angular) if(i < (this->npro *4 )) _h_ang_params[i] = dist;
        }

        //Checkers as Lambdas
        auto valid_angles = [](double a, double b, double c, double d)
        {
            double anglemin = std::min({a, b, c, d});
            double anglemax = std::max({a, b, c, d});
            if (anglemin < -1.0 || anglemax > 1.0){
                throw oxDNAException("Cos of Angle in Parameter File not in Valid bounds");
            }
        };

        auto valid_spring_params = [](int N, int x, int y, double d, char s, double k){
            if(x < 0 || x > N) throw oxDNAException("Invalid Particle ID %d in Parameter File", x);
            if(y < 0 || y > N) throw oxDNAException("Invalid Particle ID %d in Parameter File", y);
            if(d < 0) throw oxDNAException("Invalid Eq Distance %d in Parameter File", d);
            if(s != 's') throw oxDNAException("Potential Type %c Not Supported", s);
            if(k < 0) throw oxDNAException("Spring Constant %f Not Supported", k);
        };

        //Reading Parameter File
        int key1, key2 = 0;
        c_number a0, b0, c0, d0;
        std::string carbons;

        //total connections
        int spring_connection_num = 0;

        //allocate and declare affected_len vector
        _affected_len = new int[this->npro]();
        for(int i = 0; i < this->npro; i++) _affected_len[i] = 0;

        std::fstream parameters;
        parameters.open(this->_parameterfile, std::ios::in);
        getline (parameters,carbons);
        //Read Parameter File
        if (parameters.is_open())
        {
            while (parameters >> key1 >> key2 >> dist >> potswitch >> potential)
            {
                valid_spring_params(N, key1, key2, dist, potswitch, potential);
                spring_connection_num += 1;

                if(offset != 0) {
                    key1 -= offset;
                    key2 -= offset;
                }

                _affected_len[key1] += 1;
                _affected_len[key2] += 1;
                //potswitch is currently unused but may be later

                if(_angular) {
                    if (key2 - key1 == 1) {
                        //Angular Parameters
                        parameters >> a0 >> b0 >> c0 >> d0;
                        valid_angles(a0, b0, c0, d0);
                        _h_ang_params[key1 * 4] = a0;
                        _h_ang_params[key1 * 4 + 1] = b0;
                        _h_ang_params[key1 * 4 + 2] = c0;
                        _h_ang_params[key1 * 4 + 3] = d0;


                        if (this->_parameter_kbkt) {
                            parameters >> _kbend >> _ktor;
                            if (_kbend < 0 || _ktor < 0)
                                throw oxDNAException("Invalid pairwise kb/kt Value Declared in Parameter File. Check Par Formatting");
                            _h_ang_kbkt[key1 * 2] = _kbend;
                            _h_ang_kbkt[key1 * 2 + 1] = _ktor;
                        }
                    }
                }

                _spring_potential[key1*this->npro + key2] = potential;
                _spring_eqdist[key1*this->npro + key2] = dist;

                _spring_potential[key2*this->npro + key1] = potential;
                _spring_eqdist[key2*this->npro + key1] = dist;

            }
            parameters.close();
        } else {
            throw oxDNAException("ParameterFile Could Not Be Opened");
        }

        //Compressed Parameter Initialization
        _h_affected_indx = new int[this->npro + 1]();
        _h_affected = new int[spring_connection_num*2]();
        _h_aff_gamma = new c_number[spring_connection_num*2]();
        _h_aff_eqdist = new c_number[spring_connection_num*2]();
        auto zero = (c_number) 0.f;
        for(int i = 0; i < this->npro+1; i++) _h_affected_indx[i] = 0;
        for(int i = 0; i < spring_connection_num*2; i++){
            _h_affected[i] = 0;
            _h_aff_gamma[i] = zero;
            _h_aff_eqdist[i] = zero;
        }

        //Compressed Index
        int param_indx = 0;
        //For each residue
        for(int i = 0; i < this->npro; i++){
            //Fill _h_affected filtering through larger arrays filled in parameter file reading
            for(int j = i*this->npro; j < i*this->npro+this->npro; j++){
                if(_spring_eqdist[j] != 0.f){
                    //Affected List, Access is controlled with indices in _h_affected_indx
                    _h_affected[param_indx] = j % this->npro;
                    //Stored in same way for easy access, spring constants
                    _h_aff_gamma[param_indx] = _spring_potential[j];
                    //eq_distance
                    _h_aff_eqdist[param_indx] = _spring_eqdist[j];
                    param_indx += 1;
                }
            }
        }

        //Don't need Larger arrays anymore, safe to delete
        if(_spring_eqdist != NULL) delete[] _spring_eqdist;
        _spring_eqdist = NULL; //Otherwise dangling Pointer
        if(_spring_potential != NULL) delete[] _spring_potential;
        _spring_potential = NULL;

        //Allocation and Copying of Compressed Parameters
        CUDA_SAFE_CALL(cudaMalloc(&_d_affected, 2 * spring_connection_num * sizeof(int)));
        CUDA_SAFE_CALL(cudaMemcpy(_d_affected, _h_affected, 2 * spring_connection_num * sizeof(int), cudaMemcpyHostToDevice));

        CUDA_SAFE_CALL(cudaMalloc(&_d_aff_gamma, 2 * spring_connection_num * sizeof(int)));
        CUDA_SAFE_CALL(cudaMemcpy(_d_aff_gamma, _h_aff_gamma, 2 * spring_connection_num * sizeof(int), cudaMemcpyHostToDevice));

        CUDA_SAFE_CALL(cudaMalloc(&_d_aff_eqdist, 2 * spring_connection_num * sizeof(int)));
        CUDA_SAFE_CALL(cudaMemcpy(_d_aff_eqdist, _h_aff_eqdist, 2 * spring_connection_num * sizeof(int), cudaMemcpyHostToDevice));

        int ind = 0;
        _h_affected_indx[0] = 0;
        //make indx access list where: _h_affected_indx[i] lower bound of i's parameters, _h_affected_indx[i+1] upper bound of i's parameters
        for(int i = 0; i < this->npro; i++){
            ind += _affected_len[i];
            _h_affected_indx[i+1] += ind;
        }

        //Don't need this anymore
        if(_affected_len != NULL) delete[] _affected_len;
        _affected_len = NULL;

        //Allocation and copying of Indice List for accessing compressed parameters
        CUDA_SAFE_CALL(cudaMalloc(&_d_affected_indx, (this->npro+1)*sizeof(int)));
        CUDA_SAFE_CALL(cudaMemcpy(_d_affected_indx, _h_affected_indx, (this->npro+1)*sizeof(int), cudaMemcpyHostToDevice));

        if(_angular) {
            //Parameters for Bending/Torsional, _h_ang_params is filled in parameter file reading
            CUDA_SAFE_CALL(cudaMemcpy(_d_ang_params, _h_ang_params, _ang_param_size, cudaMemcpyHostToDevice));
            if (this->_parameter_kbkt) {
                CUDA_SAFE_CALL(cudaMemcpy(_d_ang_kbkt, _h_ang_kbkt, sizeof(c_number) * (this->npro * 2),
                                          cudaMemcpyHostToDevice));
            }
            float param_memory_mb;
            //Memory Used by Parameters
            if (this->_parameter_kbkt) {
                param_memory_mb =
                        (spring_connection_num * 2 * sizeof(int) + 2 * spring_connection_num * 2 * sizeof(c_number)
                         + (this->npro + 1) * sizeof(int) + 6 * this->npro * sizeof(c_number)) / SQR(1024);
            } else {
                param_memory_mb =
                        (spring_connection_num * 2 * sizeof(int) + 2 * spring_connection_num * 2 * sizeof(c_number)
                         + (this->npro + 1) * sizeof(int) + 4 * this->npro * sizeof(c_number)) / SQR(1024);
            }
            OX_LOG(Logger::LOG_INFO, "Spring Parameters Size: %.2f MB", param_memory_mb);
        }

    } else OX_LOG(Logger::LOG_INFO, "Parfile: NONE, No protein parameters were filled");

    float f_copy = this->_hb_multiplier;
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_hb_multi, &f_copy, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &N, sizeof(int)) );

    c_number tmp[50];
    for(int i = 0; i < 2; i++) for(int j = 0; j < 5; j++) for(int k = 0; k < 5; k++) tmp[i*25 + j*5 + k] = this->F1_EPS[i][j][k];

    COPY_ARRAY_TO_CONSTANT(MD_F1_EPS, tmp, 50);

    for(int i = 0; i < 2; i++) for(int j = 0; j < 5; j++) for(int k = 0; k < 5; k++) tmp[i*25 + j*5 + k] = this->F1_SHIFT[i][j][k];

    COPY_ARRAY_TO_CONSTANT(MD_F1_SHIFT, tmp, 50);

    COPY_ARRAY_TO_CONSTANT(MD_F1_A, this->F1_A, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F1_RC, this->F1_RC, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F1_R0, this->F1_R0, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F1_BLOW, this->F1_BLOW, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F1_BHIGH, this->F1_BHIGH, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F1_RLOW, this->F1_RLOW, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F1_RHIGH, this->F1_RHIGH, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F1_RCLOW, this->F1_RCLOW, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F1_RCHIGH, this->F1_RCHIGH, 2);

    COPY_ARRAY_TO_CONSTANT(MD_F2_K, this->F2_K, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F2_RC, this->F2_RC, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F2_R0, this->F2_R0, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F2_BLOW, this->F2_BLOW, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F2_BHIGH, this->F2_BHIGH, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F2_RLOW, this->F2_RLOW, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F2_RHIGH, this->F2_RHIGH, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F2_RCLOW, this->F2_RCLOW, 2);
    COPY_ARRAY_TO_CONSTANT(MD_F2_RCHIGH, this->F2_RCHIGH, 2);

    COPY_ARRAY_TO_CONSTANT(MD_F5_PHI_A, this->F5_PHI_A, 4);
    COPY_ARRAY_TO_CONSTANT(MD_F5_PHI_B, this->F5_PHI_B, 4);
    COPY_ARRAY_TO_CONSTANT(MD_F5_PHI_XC, this->F5_PHI_XC, 4);
    COPY_ARRAY_TO_CONSTANT(MD_F5_PHI_XS, this->F5_PHI_XS, 4);


    if(this->_use_edge) CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_n_forces, &this->_n_forces, sizeof(int)));

    // copied from DNA2Interaction::init() (CPU), the least bad way of doing things
    // We wish to normalise with respect to T=300K, I=1M. 300K=0.1 s.u. so divide this->_T by 0.1
    c_number lambda = _debye_huckel_lambdafactor * sqrt(this->_T / 0.1f) / sqrt(_salt_concentration);
    // RHIGH gives the distance at which the smoothing begins
    _debye_huckel_RHIGH = 3.0 * lambda;
    _minus_kappa = -1.0 / lambda;

    // these are just for convenience for the smoothing parameter computation
    c_number x = _debye_huckel_RHIGH;
    c_number q = _debye_huckel_prefactor;
    c_number l = lambda;

    // compute the some smoothing parameters
    _debye_huckel_B = -(exp(-x / l) * q * q * (x + l) * (x + l)) / (-4. * x * x * x * l * l * q);
    _debye_huckel_RC = x * (q * x + 3. * q * l) / (q * (x + l));

    c_number debyecut;
    if (this->_grooving) {
        debyecut = 2.0f * sqrt(SQR(POS_MM_BACK1) + SQR(POS_MM_BACK2)) + _debye_huckel_RC;
    } else {
        debyecut = 2.0f * sqrt(SQR(POS_BACK)) + _debye_huckel_RC;
    }
    // the cutoff radius for the potential should be the larger of rcut and debyecut
    if (debyecut > this->_rcut) {
        this->_rcut = debyecut;
        this->_sqr_rcut = debyecut * debyecut;
    }
    // End copy from DNA2Interaction

    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_dh_RC, &_debye_huckel_RC, sizeof(float)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_dh_RHIGH, &_debye_huckel_RHIGH, sizeof(float)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_dh_prefactor, &_debye_huckel_prefactor, sizeof(float)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_dh_B, &_debye_huckel_B, sizeof(float)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_dh_minus_kappa, &_minus_kappa, sizeof(float)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_dh_half_charged_ends, &_debye_huckel_half_charged_ends, sizeof(bool)));

    //Constants for DNA/Protein Excluded Volume Interactions
    //Backbone-Protein Excluded Volume Parameters
    _pro_backbone_sigma = 0.57f;
    _pro_backbone_rstar= 0.569f;
    _pro_backbone_b = 178699253.5f;
    _pro_backbone_rcut = 0.572934f;
    _pro_backbone_sqr_rcut = 0.3283f;
    //Base-Protein Excluded Volume Parameters
    _pro_base_sigma = 0.36f;
    _pro_base_rstar= 0.359f;
    _pro_base_b = 296866090.0f;
    _pro_base_rcut = 0.362897f;
    _pro_base_sqr_rcut = 0.1317f;
    //Protein-Protein Excluded Volume Parameters
    _pro_sigma = 0.35f;
    _pro_rstar = 0.349f;
    _pro_b = 306484596.0f;
    _pro_rcut = 0.352894;
    _pro_sqr_rcut = 0.12454f; //_rc ^2

    if(_angular) {
        if (!this->_parameter_kbkt) {
            _kbend = this->_k_bend;
            _ktor = this->_k_tor;

            //kb and kt Parameters
            CUDA_SAFE_CALL(cudaMemcpyToSymbol(_kb, &_kbend, sizeof(float)));
            CUDA_SAFE_CALL(cudaMemcpyToSymbol(_kt, &_ktor, sizeof(float)));
        }
    }

    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_sigma, &_pro_sigma, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_rstar, &_pro_rstar, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_rc, &_pro_rcut, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_b, &_pro_b, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_sqr_rcut, &_pro_sqr_rcut, sizeof(float)) );

    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_backbone_sigma, &_pro_backbone_sigma, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_backbone_rstar, &_pro_backbone_rstar, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_backbone_rc, &_pro_backbone_rcut, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_backbone_b, &_pro_backbone_b, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_backbone_sqr_rcut, &_pro_backbone_sqr_rcut, sizeof(float)) );

    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_base_sigma, &_pro_base_sigma, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_base_rstar, &_pro_base_rstar, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_base_rc, &_pro_base_rcut, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_base_b, &_pro_base_b, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_pro_base_sqr_rcut, &_pro_base_sqr_rcut, sizeof(float)) );

    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dna_sqr_rcut, &this->_sqr_rcut, sizeof(float)) );

    //Parameters for DNANM book keeping
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(_ndna, &this->ndna, sizeof(int)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(_npro, &this->npro, sizeof(int)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(_offset, &this->offset, sizeof(int)) );
}

void CUDADNANMInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
    if(_d_is_strand_end == nullptr) {
        _init_strand_ends(d_bonds);
    }


    if(_update_st) {
        CUDA_SAFE_CALL(cudaMemset(_d_st, 0, _N * sizeof(CUDAStressTensor)));
    }

    if(_use_edge) {

        if(_n_forces == 1) { // we can directly use d_forces and d_torques so that no sum is required
            dnanm_forces_edge_nonbonded
            <<<(lists->N_edges - 1)/(_launch_cfg.threads_per_block) + 1, _launch_cfg.threads_per_block>>>
                    (d_poss, d_orientations, d_forces, d_torques, lists->d_edge_list, lists->N_edges, _d_is_strand_end,
                     _grooving, _use_debye_huckel, _use_oxDNA2_coaxial_stacking, _update_st, _d_st, d_box);
        }
        else { // sum required, somewhat slower
            dnanm_forces_edge_nonbonded
            <<<(lists->N_edges - 1)/(_launch_cfg.threads_per_block) + 1, _launch_cfg.threads_per_block>>>
                    (d_poss, d_orientations, _d_edge_forces, _d_edge_torques, lists->d_edge_list, lists->N_edges, _d_is_strand_end,
                     _grooving, _use_debye_huckel, _use_oxDNA2_coaxial_stacking, _update_st, _d_st, d_box);

            _sum_edge_forces_torques(d_forces, d_torques);
        }

        dna_forces_edge_bonded_dnanm
        <<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
                (d_poss, d_orientations, d_forces, d_torques, d_bonds, _grooving, _use_oxDNA2_FENE, _use_mbf, _mbf_xmax,
                 _mbf_finf, _update_st, _d_st);

        protein_forces_edge_bonded
        <<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
                (d_poss, d_orientations, d_forces, d_torques, d_bonds, _grooving, _use_oxDNA2_FENE, _use_mbf, _mbf_xmax,
                 _mbf_finf, d_box, _update_st, _d_st, _d_aff_eqdist, _d_aff_gamma, _d_affected_indx, _d_affected);
    } else throw oxDNAException("Edge Approach is only implemented for DNANM Interaction using CUDA approach. Please add use_edge = 1 to your input file.");


//            if(_angular){
//                dnanm_forces_edge_bonded_angular
//                        <<< this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block >>>
//                (d_poss, d_orientations, d_forces, d_torques, d_bonds, this->_grooving, _use_oxDNA2_FENE, this->_use_mbf, this->_mbf_xmax, this->_mbf_finf, d_box, _d_aff_eqdist, _d_aff_gamma, _d_ang_params, _d_ang_kbkt, _d_affected_indx, _d_affected);
//
//            } else {
//                dnanm_forces_edge_bonded
//                        <<< this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block >>>
//                (d_poss, d_orientations, d_forces, d_torques, d_bonds, this->_grooving, _use_oxDNA2_FENE, this->_use_mbf, this->_mbf_xmax, this->_mbf_finf, d_box, _d_aff_eqdist, _d_aff_gamma, _d_affected_indx, _d_affected);
//            }

//    } else throw oxDNAException("Must Use with Lists to run simulation");
}

void CUDADNANMInteraction::_on_T_update() {
    cuda_init(_N);
}

void CUDADNANMInteraction::_init_strand_ends(LR_bonds *d_bonds) {
    CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&_d_is_strand_end, sizeof(int) * _N));
    dna_init_DNA_strand_ends<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>(_d_is_strand_end, d_bonds, _N);
}
