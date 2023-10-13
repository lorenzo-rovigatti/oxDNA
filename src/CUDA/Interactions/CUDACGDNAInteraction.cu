//
// Created by jonah on 3/26/21.
//


#include "CUDACGDNAInteraction.h"
#include "CUDA_CGDNA.cuh"
#include "../Lists/CUDASimpleVerletList.h"
#include "../Lists/CUDANoList.h"
#include "../../Interactions/CGDNAInteraction.h"

#include <sstream>
#include <fstream>


CUDACGDNAInteraction::CUDACGDNAInteraction(bool btp) : CGDNAInteraction(btp), CUDABaseInteraction() {
    
    _read_par = false;
    _angular = btp;
    //Not copied over to device memory
    _gs_spring_potential = NULL;
    _gs_spring_eqdist = NULL;
    _gs_affected_len = NULL;

    _pro_spring_potential = NULL;
    _pro_spring_eqdist = NULL;
    _pro_affected_len = NULL;

    //Copied over to device memory
    _h_pro_affected_indx = NULL;
    _d_pro_affected_indx = NULL;

    _h_gs_affected_indx = NULL;
    _d_gs_affected_indx = NULL;

    _h_pro_affected = NULL;
    _d_pro_affected = NULL;

    _h_gs_affected = NULL;
    _d_gs_affected = NULL;

    _h_pro_aff_eqdist = NULL;
    _d_pro_aff_eqdist = NULL;

    _h_gs_aff_eqdist = NULL;
    _d_gs_aff_eqdist = NULL;

    _h_pro_aff_gamma = NULL;
    _d_pro_aff_gamma = NULL;

    _h_gs_aff_gamma = NULL;
    _d_gs_aff_gamma = NULL;

    _pro_spring_param_size_number = 0;
    _gs_spring_param_size_number = 0;

    _gs_offset = -1; // set to 0 or greater if gs in system
    _pro_offset = -1; // set to 0 or greater if pro in system

    _h_gs_gs_exc_vol_params = NULL;
    _d_gs_gs_exc_vol_params = NULL;

    _h_gs_other_exc_vol_params = NULL;
    _d_gs_other_exc_vol_params = NULL;

    _edge_compatible = true;
}


CUDACGDNAInteraction::~CUDACGDNAInteraction() {
    //Delete All pointers required for spring potential and other parameters
    if(_pro_spring_potential != NULL) delete[] _pro_spring_potential;
    if(_pro_spring_eqdist != NULL) delete[] _pro_spring_eqdist;

    if(_gs_spring_potential != NULL) delete[] _gs_spring_potential;
    if(_gs_spring_eqdist != NULL) delete[] _gs_spring_eqdist;

    if(_pro_affected_len != NULL) delete[] _pro_affected_len;
    if(_gs_affected_len != NULL) delete[] _gs_affected_len;

    if(_d_pro_affected != NULL) CUDA_SAFE_CALL(cudaFree(_d_pro_affected));
    if(_d_pro_aff_gamma != NULL) CUDA_SAFE_CALL(cudaFree(_d_pro_aff_gamma));
    if(_d_pro_aff_eqdist != NULL) CUDA_SAFE_CALL(cudaFree(_d_pro_aff_eqdist));
    if(_d_pro_affected_indx != NULL) CUDA_SAFE_CALL( cudaFree(_d_pro_affected_indx) );

    if(_d_gs_affected != NULL) CUDA_SAFE_CALL(cudaFree(_d_gs_affected));
    if(_d_gs_aff_gamma != NULL) CUDA_SAFE_CALL(cudaFree(_d_gs_aff_gamma));
    if(_d_gs_aff_eqdist != NULL) CUDA_SAFE_CALL(cudaFree(_d_gs_aff_eqdist));
    if(_d_gs_affected_indx != NULL) CUDA_SAFE_CALL( cudaFree(_d_gs_affected_indx) );

    if(_d_gs_other_exc_vol_params != NULL) CUDA_SAFE_CALL( cudaFree(_d_gs_other_exc_vol_params));
    if(_d_gs_gs_exc_vol_params != NULL) CUDA_SAFE_CALL( cudaFree(_d_gs_gs_exc_vol_params));

    if(_h_pro_affected != NULL) delete[] _h_pro_affected;
    if(_h_pro_aff_gamma != NULL) delete[] _h_pro_aff_gamma;
    if(_h_pro_aff_eqdist != NULL) delete[] _h_pro_aff_eqdist;
    if(_h_pro_affected_indx != NULL) delete[] _h_pro_affected_indx;

    if(_h_gs_affected != NULL) delete[] _h_gs_affected;
    if(_h_gs_aff_gamma != NULL) delete[] _h_gs_aff_gamma;
    if(_h_gs_aff_eqdist != NULL) delete[] _h_gs_aff_eqdist;
    if(_h_gs_affected_indx != NULL) delete[] _h_gs_affected_indx;

    if(_h_gs_other_exc_vol_params != NULL) delete[] _h_gs_other_exc_vol_params;
    if(_h_gs_gs_exc_vol_params != NULL) delete[] _h_gs_gs_exc_vol_params;
}


void CUDACGDNAInteraction::get_settings(input_file &inp) {
    std::string inter_type;
    if (getInputString(&inp, "parfile", _parameterfile, 0) != KEY_FOUND) {
        throw oxDNAException("Key 'parfile' not found. Necessary for Protein sims.");
    }

    char s[5] = "none";
    if(strcmp(this->_parameterfile, s) != 0) _read_par = true;

    if (!getInputString(&inp, "topology", this->_topology_filename, 0) == KEY_FOUND){
        throw oxDNAException("Key 'topology_file' not found.");
    }
    CGDNAInteraction::get_settings(inp);
}


void CUDACGDNAInteraction::cuda_init(int N) {
    CUDABaseInteraction::cuda_init(N);
    CGDNAInteraction::init();

    std::fstream top;
    int tmp1, tmp2;
    char line[5120];

    std::string base_tmp;
    std::array<bool, 3> found = {false, false, false};
    std::array<int, 3> element_numbers = {0, 0, 0};

    int gs_species = 0;
    int subtype, strand;

    top.open(this->_topology_filename, std::ios::in);
    if (top.is_open()) {
        top >> tmp1 >> tmp2 >> this->ndna >> this->npro >> this->ndnas >> this->npep;
        this->ngstrands = tmp2 - this->ndnas - this->npep;
        this->ngs = tmp1 - this->ndna - this->npro;  // number of gs particles

        while (top.good()) {
            top.getline(line, 5120);
            if (strlen(line) == 0 || line[0] == '#')
                continue;

            std::stringstream ss(line);
            ss >> strand >> base_tmp;
            // CG DNA flag
            // negative strand number and contains digits in the base field "gs0"
            if (strand < 0 && !found[0] && std::any_of(base_tmp.begin(), base_tmp.end(), ::isdigit)) {
                this->topology_order[std::count(found.begin(), found.end(), true)] = 0;
                element_numbers[std::count(found.begin(), found.end(), true)] = this->ngs;
                found[0] = true;

            }
            // Protein Flag
            // negative strand number and no digits in the base field -1 T
            if (strand < 0 && !found[1] && !std::any_of(base_tmp.begin(), base_tmp.end(), ::isdigit)) {
                this->topology_order[std::count(found.begin(), found.end(), true)] = 1;
                element_numbers[std::count(found.begin(), found.end(), true)] = this->npro;
                found[1] = true;
            }
            // DNA Flag
            if (strand > 0 && !found[2]) { // positive strand number DNA 1 A
                this->topology_order[std::count(found.begin(), found.end(), true)] = 2;
                element_numbers[std::count(found.begin(), found.end(), true)] = this->ndna;
                found[2] = true;
            }


            // need to find the number of gs particle types for parameters (system specific)
            if(strand < 0){
                if(std::any_of(base_tmp.begin(), base_tmp.end(),::isdigit)){ // if digit in base field it's a gs particle
                    base_tmp.erase(0, 2); // remove first 2 characters -> 'gs', left with subtype integer
                    subtype=stoi(base_tmp);
                    if(subtype> gs_species){
                        gs_species = subtype;
                    }
                }
            }

        }

    } else {
        throw oxDNAException("Could not open Topology File");
    }

    int current_offset = 0;
    for(auto i = 0; i < std::count(found.begin(), found.end(), true); i++){
        if(topology_order[i] == 0){ // gs particle
            _gs_offset = current_offset;
        } else if(topology_order[i] == 1) {
            _pro_offset = current_offset;
        }

        current_offset += element_numbers[i];
    }

    if(_read_par){
        //Initalizing Some Host and Device Arrays for Spring Parameters
        _pro_spring_param_size_number = sizeof(c_number) * (this->npro*this->npro);
        _gs_spring_param_size_number = sizeof(c_number) * (this->ngs*this->ngs);
        //_ang_param_size = sizeof(c_number) * (this->npro*4);

        _pro_spring_potential = new c_number[this->npro*this->npro]();
        _pro_spring_eqdist = new c_number[this->npro*this->npro]();

        _gs_spring_potential = new c_number[this->ngs*this->ngs]();
        _gs_spring_eqdist = new c_number[this->ngs*this->ngs]();

        char potswitch = 'x';
        c_number potential = 0.f;
        c_number dist = 0.f;
        //initializing array members
        for(int i = 0; i< (this->npro*this->npro); i++){
            _pro_spring_eqdist[i] = dist;
            _pro_spring_potential[i] = potential;
        }
        for(int i = 0; i< (this->ngs*this->ngs); i++){
            _gs_spring_eqdist[i] = dist;
            _gs_spring_potential[i] = potential;
        }

        auto valid_spring_params = [](int N, int x, int y, double d, char s, double k){
            if(x < 0 || x > N) throw oxDNAException("Invalid Particle ID %d in Parameter File", x);
            if(y < 0 || y > N) throw oxDNAException("Invalid Particle ID %d in Parameter File", y);
            if(d < 0) throw oxDNAException("Invalid Eq Distance %d in Parameter File", d);
            if(s != 's') throw oxDNAException("Potential Type %c Not Supported", s);
            if(k < 0) throw oxDNAException("Spring Constant %f Not Supported", k);
        };

        //Reading Parameter File
        int key1, key2 = 0;
        //c_number a0, b0, c0, d0;
        std::string carbons;

        //total connections
        int pro_spring_connection_num = 0;
        int gs_spring_connection_num = 0;

        //allocate and declare affected_len vector
        _pro_affected_len = new int[this->npro]();
        for(int i = 0; i < this->npro; i++) _pro_affected_len[i] = 0;

        _gs_affected_len = new int[this->ngs]();
        for(int i = 0; i < this->ngs; i++) _gs_affected_len[i] = 0;

        std::fstream parameters;
        parameters.open(this->_parameterfile, std::ios::in);
        getline (parameters,carbons);
        //Read Parameter File
        if (parameters.is_open())
        {
            while (parameters >> key1 >> key2 >> dist >> potswitch >> potential)
            {
                valid_spring_params(N, key1, key2, dist, potswitch, potential);


                if(key1 < _gs_offset || _gs_offset == -1){ // has to be a protein spring parameter
                    key1 -= _pro_offset;
                    key2 -= _pro_offset;

                    pro_spring_connection_num += 1;

                    _pro_affected_len[key1] += 1;
                    _pro_affected_len[key2] += 1;

                    _pro_spring_potential[key1*this->npro + key2] = potential;
                    _pro_spring_eqdist[key1*this->npro + key2] = dist;

                    _pro_spring_potential[key2*this->npro + key1] = potential;
                    _pro_spring_eqdist[key2*this->npro + key1] = dist;

                } else if (key1 < _pro_offset || _pro_offset == -1){ // has to be a gs spring parameter
                    key1 -= _gs_offset;
                    key2 -= _gs_offset;

                    gs_spring_connection_num += 1;

                    _gs_affected_len[key1] += 1;
                    _gs_affected_len[key2] += 1;

                    _gs_spring_potential[key1*this->ngs + key2] = potential;
                    _gs_spring_eqdist[key1*this->ngs + key2] = dist;

                    _gs_spring_potential[key2*this->ngs + key1] = potential;
                    _gs_spring_eqdist[key2*this->ngs + key1] = dist;

                } else {
                    throw oxDNAException("Particle Id %f Not Compatible with system", key1);
                }

            }
            parameters.close();
        } else {
            throw oxDNAException("ParameterFile Could Not Be Opened");
        }

        //Compressed Parameter Initialization
        _h_pro_affected_indx = new int[this->npro + 1]();
        _h_pro_affected = new int[pro_spring_connection_num*2]();
        _h_pro_aff_gamma = new c_number[pro_spring_connection_num*2]();
        _h_pro_aff_eqdist = new c_number[pro_spring_connection_num*2]();
        auto zero = (c_number) 0.f;
        for(int i = 0; i < this->npro+1; i++) _h_pro_affected_indx[i] = 0;
        for(int i = 0; i < pro_spring_connection_num*2; i++){
            _h_pro_affected[i] = 0;
            _h_pro_aff_gamma[i] = zero;
            _h_pro_aff_eqdist[i] = zero;
        }

        //Compressed Index
        int param_indx = 0;
        //For each residue
        for(int i = 0; i < this->npro; i++){
            //Fill _h_affected filtering through larger arrays filled in parameter file reading
            for(int j = i*this->npro; j < i*this->npro+this->npro; j++){
                if(_pro_spring_eqdist[j] != 0.f){
                    //Affected List, Access is controlled with indices in _h_affected_indx
                    _h_pro_affected[param_indx] = j % this->npro;
                    //Stored in same way for easy access, spring constants
                    _h_pro_aff_gamma[param_indx] = _pro_spring_potential[j];
                    //eq_distance
                    _h_pro_aff_eqdist[param_indx] = _pro_spring_eqdist[j];
                    param_indx += 1;
                }
            }
        }

        //Don't need Larger arrays anymore, safe to delete
        if(_pro_spring_eqdist != NULL) delete[] _pro_spring_eqdist;
        _pro_spring_eqdist = NULL; //Otherwise dangling Pointer
        if(_pro_spring_potential != NULL) delete[] _pro_spring_potential;
        _pro_spring_potential = NULL;

        //Allocation and Copying of Compressed Parameters
        CUDA_SAFE_CALL(cudaMalloc(&_d_pro_affected, 2 * pro_spring_connection_num * sizeof(int)));
        CUDA_SAFE_CALL(cudaMemcpy(_d_pro_affected, _h_pro_affected, 2 * pro_spring_connection_num * sizeof(int), cudaMemcpyHostToDevice));

        CUDA_SAFE_CALL(cudaMalloc(&_d_pro_aff_gamma, 2 * pro_spring_connection_num * sizeof(c_number)));
        CUDA_SAFE_CALL(cudaMemcpy(_d_pro_aff_gamma, _h_pro_aff_gamma, 2 * pro_spring_connection_num * sizeof(c_number), cudaMemcpyHostToDevice));

        CUDA_SAFE_CALL(cudaMalloc(&_d_pro_aff_eqdist, 2 * pro_spring_connection_num * sizeof(c_number)));
        CUDA_SAFE_CALL(cudaMemcpy(_d_pro_aff_eqdist, _h_pro_aff_eqdist, 2 * pro_spring_connection_num * sizeof(c_number), cudaMemcpyHostToDevice));

        int ind = 0;
        _h_pro_affected_indx[0] = 0;
        //make indx access list where: _h_affected_indx[i] lower bound of i's parameters, _h_affected_indx[i+1] upper bound of i's parameters
        for(int i = 0; i < this->npro; i++){
            ind += _pro_affected_len[i];
            _h_pro_affected_indx[i+1] += ind;
        }

        //Don't need this anymore
        if(_pro_affected_len != NULL) delete[] _pro_affected_len;
        _pro_affected_len = NULL;

        //Allocation and copying of Indice List for accessing compressed parameters
        CUDA_SAFE_CALL(cudaMalloc(&_d_pro_affected_indx, (this->npro+1)*sizeof(int)));
        CUDA_SAFE_CALL( cudaMemcpy(_d_pro_affected_indx, _h_pro_affected_indx, (this->npro+1)*sizeof(int), cudaMemcpyHostToDevice));

        /* ///////////////////////////////////////////////// */

        // Same exact thing as above, but for gs
        //Compressed Parameter Initialization
        _h_gs_affected_indx = new int[this->ngs + 1]();
        _h_gs_affected = new int[gs_spring_connection_num*2]();
        _h_gs_aff_gamma = new c_number[gs_spring_connection_num*2]();
        _h_gs_aff_eqdist = new c_number[gs_spring_connection_num*2]();
        for(int i = 0; i < this->ngs+1; i++) _h_gs_affected_indx[i] = 0;
        for(int i = 0; i < gs_spring_connection_num*2; i++){
            _h_gs_affected[i] = 0;
            _h_gs_aff_gamma[i] = zero;
            _h_gs_aff_eqdist[i] = zero;
        }

        //Compressed Index
        param_indx = 0;
        //For each residue
        for(int i = 0; i < this->ngs; i++){
            //Fill _h_affected filtering through larger arrays filled in parameter file reading
            for(int j = i*this->ngs; j < i*this->ngs+this->ngs; j++){
                if(_gs_spring_eqdist[j] != 0.f){
                    //Affected List, Access is controlled with indices in _h_affected_indx
                    _h_gs_affected[param_indx] = j % this->ngs;
                    //Stored in same way for easy access, spring constants
                    _h_gs_aff_gamma[param_indx] = _gs_spring_potential[j];
                    //eq_distance
                    _h_gs_aff_eqdist[param_indx] = _gs_spring_eqdist[j];
                    param_indx += 1;
                }
            }
        }

        //Don't need Larger arrays anymore, safe to delete
        if(_gs_spring_eqdist != NULL) delete[] _gs_spring_eqdist;
        _gs_spring_eqdist = NULL; //Otherwise dangling Pointer
        if(_gs_spring_potential != NULL) delete[] _gs_spring_potential;
        _gs_spring_potential = NULL;

        //Allocation and Copying of Compressed Parameters
        CUDA_SAFE_CALL(cudaMalloc(&_d_gs_affected, 2 * gs_spring_connection_num * sizeof(int)));
        CUDA_SAFE_CALL(cudaMemcpy(_d_gs_affected, _h_gs_affected, 2 * gs_spring_connection_num * sizeof(int), cudaMemcpyHostToDevice));

        CUDA_SAFE_CALL(cudaMalloc(&_d_gs_aff_gamma, 2 * gs_spring_connection_num * sizeof(c_number)));
        CUDA_SAFE_CALL(cudaMemcpy(_d_gs_aff_gamma, _h_gs_aff_gamma, 2 * gs_spring_connection_num * sizeof(c_number), cudaMemcpyHostToDevice));

        CUDA_SAFE_CALL(cudaMalloc(&_d_gs_aff_eqdist, 2 * gs_spring_connection_num * sizeof(c_number)));
        CUDA_SAFE_CALL(cudaMemcpy(_d_gs_aff_eqdist, _h_gs_aff_eqdist, 2 * gs_spring_connection_num * sizeof(c_number), cudaMemcpyHostToDevice));

        ind = 0;
        _h_gs_affected_indx[0] = 0;
        //make indx access list where: _h_affected_indx[i] lower bound of i's parameters, _h_affected_indx[i+1] upper bound of i's parameters
        for(int i = 0; i < this->ngs; i++){
            ind += _gs_affected_len[i];
            _h_gs_affected_indx[i+1] += ind;
        }

        //Don't need this anymore
        if(_gs_affected_len != NULL) delete[] _gs_affected_len;
        _gs_affected_len = NULL;

        //Allocation and copying of Indice List for accessing compressed parameters
        CUDA_SAFE_CALL(cudaMalloc(&_d_gs_affected_indx, (this->ngs+1)*sizeof(int)));
        CUDA_SAFE_CALL( cudaMemcpy(_d_gs_affected_indx, _h_gs_affected_indx, (this->ngs+1)*sizeof(int), cudaMemcpyHostToDevice));



        float param_memory_mb;
        int spring_connection_num = gs_spring_connection_num + pro_spring_connection_num;

        // eqdist and k for each spring
        param_memory_mb = (spring_connection_num * 2 * sizeof(int) + 2 * spring_connection_num * 2 * sizeof(c_number)
                 + ((this->npro + this->ngs) + 1) * sizeof(int)) / SQR(1024);

        OX_LOG(Logger::LOG_INFO, "Spring Parameters Size: %.2f MB", param_memory_mb);

    } else OX_LOG(Logger::LOG_INFO, "Parfile: NONE, No protein parameters were filled");

    /* DNA Parameters Now                          */

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


    if(this->_use_edge) CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_n_forces, &this->_n_forces, sizeof(int)) );

    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_RC, &_debye_huckel_RC, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_RHIGH, &_debye_huckel_RHIGH, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_prefactor, &_debye_huckel_prefactor, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_B, &_debye_huckel_B, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_minus_kappa, &_minus_kappa, sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dh_half_charged_ends, &_debye_huckel_half_charged_ends, sizeof(bool)) );

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

    //Parameters for CGDNA book keeping
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(_ndna, &this->ndna, sizeof(int)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(_npro, &this->npro, sizeof(int)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(_ngs, &this->ngs, sizeof(int)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(pro_offset, &_pro_offset, sizeof(int)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(gs_offset, &_gs_offset, sizeof(int)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(_gs_species, &gs_species, sizeof(int)) );


    /* Now is the fun part, making all the excluded volume parameters! Specific to GS */

    int gs_exc_vol_size = SQR(gs_species);
    int G = gs_species;
    _h_gs_gs_exc_vol_params = new c_number[gs_exc_vol_size * 4]();
    for(int i = 0; i<gs_exc_vol_size * 4; i++) _h_gs_gs_exc_vol_params[i] = 0.f;

    number max_rcut = 0.f; // max distance b/t nonbonded particles to consider
    number sigma, rstar, rc, b, interaction_dist;
    bool success;

    // fill in parameters for each pair of radii
    for (int i = 0; i < G; i++) {
        for (int j = 0; j < G; j++) {
            if (i > j)
                continue;
            interaction_dist = this->radii[27 + i] + this->radii[27 + j];
            success = _fill_in_constants(interaction_dist, sigma, rstar, b, rc);
            if(!success) throw oxDNAException("No smooth parameters found for gs particles gs types %d and %d", i, j);

            // for debug
            //printf("interaction dist %.3f \n", interaction_dist);
            //printf("sigma %.3f rstar %.3f b %.3f rc %.3f \n", sigma, rstar, b, rc);

            if(rc > max_rcut) max_rcut = rc;
            // i*gs_subtype_num-(i-1)+j
            _h_gs_gs_exc_vol_params[4 * (i*G + j)] = sigma;
            _h_gs_gs_exc_vol_params[4 * (i*G + j) + 1] = rstar;
            _h_gs_gs_exc_vol_params[4 * (i*G + j) + 2] = b;
            _h_gs_gs_exc_vol_params[4 * (i*G + j) + 3] = rc;

            _h_gs_gs_exc_vol_params[4 * (j*G + i)] = sigma;
            _h_gs_gs_exc_vol_params[4 * (j*G + i) + 1] = rstar;
            _h_gs_gs_exc_vol_params[4 * (j*G + i) + 2] = b;
            _h_gs_gs_exc_vol_params[4 * (j*G + i) + 3] = rc;
        }
    }

    int excl_vol_pro_dna_size = G * 3;
    _h_gs_other_exc_vol_params = new c_number[excl_vol_pro_dna_size * 4]();
    for(int i = 0; i<excl_vol_pro_dna_size * 4; i++) _h_gs_other_exc_vol_params[i] = 0.f;
    // fill in protein and dna excl volume
    for (int i = 0; i < G; i++) {
        for (int j = 0; j < 3; j++) {
            if (j == 0) {
                _fill_in_constants(radii[i+27] + _pro_sigma / 2, sigma, rstar, b, rc); // protein 0
            } else if (j == 1) {
                _fill_in_constants(radii[i+27] + _pro_base_sigma - .175f, sigma, rstar, b, rc); // dna base 1
            } else if (j == 2) {
                _fill_in_constants(radii[i+27] + _pro_backbone_sigma - .175f, sigma, rstar, b, rc); // dna backbone 2
            } else {
                throw oxDNAException("Problem filling in Excluded Volume Constants");
            }

            if(rc > max_rcut) max_rcut = rc;

            _h_gs_other_exc_vol_params[4 * (3 * i + j)] = sigma;
            _h_gs_other_exc_vol_params[4 * (3 * i + j) + 1] = rstar;
            _h_gs_other_exc_vol_params[4 * (3 * i + j) + 2] = b;
            _h_gs_other_exc_vol_params[4 * (3 * i + j) + 3] = rc;
        }
    }

    CUDA_SAFE_CALL(cudaMalloc(&_d_gs_gs_exc_vol_params, gs_exc_vol_size * 4 * sizeof(c_number)));
    CUDA_SAFE_CALL(cudaMemcpy(_d_gs_gs_exc_vol_params, _h_gs_gs_exc_vol_params, gs_exc_vol_size * 4 * sizeof(c_number), cudaMemcpyHostToDevice));

    CUDA_SAFE_CALL(cudaMalloc(&_d_gs_other_exc_vol_params, excl_vol_pro_dna_size * 4 * sizeof(c_number)));
    CUDA_SAFE_CALL(cudaMemcpy(_d_gs_other_exc_vol_params, _h_gs_other_exc_vol_params, excl_vol_pro_dna_size * 4 * sizeof(c_number), cudaMemcpyHostToDevice));

    _rcut = max_rcut; // set the rcut. Used by backend verlet lists
    /* Now all Parameters are copied to the device so we can go ahead and run the model */
}

void CUDACGDNAInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
    CUDASimpleVerletList *_v_lists = dynamic_cast<CUDASimpleVerletList *>(lists);
    if(_v_lists != NULL) {
        if (_v_lists->use_edge()) {

            cgdna_forces_edge_nonbonded
            <<<(_v_lists->N_edges - 1) / (this->_launch_cfg.threads_per_block) +
               1, this->_launch_cfg.threads_per_block>>>
                    (d_poss, d_orientations, this->_d_edge_forces, this->_d_edge_torques, _v_lists->d_edge_list,
                     _v_lists->N_edges, d_bonds, this->_grooving, _use_debye_huckel, _use_oxDNA2_coaxial_stacking,
                     d_box, _d_gs_gs_exc_vol_params, _d_gs_other_exc_vol_params);

            this->_sum_edge_forces_torques(d_forces, d_torques);

            // potential for removal here
            cudaThreadSynchronize();
            CUT_CHECK_ERROR("forces_second_step error -- after non-bonded");


            cgdna_forces_edge_bonded
            <<< this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block >>>
            (d_poss, d_orientations, d_forces, d_torques, d_bonds, this->_grooving, _use_oxDNA2_FENE, this->_use_mbf,
             this->_mbf_xmax, this->_mbf_finf, d_box, _d_pro_aff_eqdist, _d_gs_aff_eqdist, _d_pro_aff_gamma, _d_gs_aff_gamma,
             _d_pro_affected_indx, _d_gs_affected_indx, _d_pro_affected, _d_gs_affected);


        } else throw oxDNAException("Edge Approach is only implemented for DNANM Interaction using CUDA approach. Please add use_edge = 1 to your input file.");

    } else throw oxDNAException("Must Use with Lists to run simulation");
}

void CUDACGDNAInteraction::_on_T_update() {
    cuda_init(_N);
}
