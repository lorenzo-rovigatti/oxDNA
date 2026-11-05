//
// Created by jonah on 3/26/21.
//


#include "CUDADNANMInteraction.h"
#include "CUDA_DNANM.cuh"
#include "../../Interactions/DNANMInteraction.h"

#include <algorithm>
#include <fstream>

CUDADNANMInteraction::CUDADNANMInteraction(bool btp) : DNANMInteraction(btp), CUDABaseInteraction() {
    
    _read_par = false;
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
    // the bending/torsion kernel (dnanm_forces_edge_bonded_angular) has not been ported to the edge-based force
    // computation, so DNACT would run with DNANM dynamics while reporting DNACT energies
    if(this->_anm.torsion()) {
        throw oxDNAException("DNACT is not supported on CUDA yet (the bending and torsion forces are not implemented on the GPU): use backend = CPU, or interaction_type = DNANM");
    }

    DNANMInteraction::get_settings(inp);
    _read_par = this->_anm.has_parameter_file();

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
//    this object is not the one whose read_topology is called by the backend, so ndna, npro, ndnas and _firststrand
//    have to be read here (old- or new-style topology)
    this->read_topology_counts();
    _n_nucleotides = this->_anm_topology.N_nucleotides;
    _n_amino_acids = this->_anm_topology.N_amino_acids;
    // the spring arrays are indexed by the index of the amino acid within the protein block
    offset = this->_anm_topology.protein_first ? 0 : _n_nucleotides;


    if(_read_par) {
        // the parameter file is parsed by the CPU protein model; here the springs are stored in a compressed format:
        // the springs of the amino acid with local index i (i.e. particle index - offset) are those in
        // [_h_affected_indx[i], _h_affected_indx[i + 1]), sorted by the local index of the other amino acid
        this->_anm.parse_parameter_file(N);
        const auto &springs = this->_anm.springs();
        int n_pro = _n_amino_acids;

        std::vector<std::vector<std::pair<int, const ANMModel::Spring *>>> per_residue(n_pro);
        for(auto &spring : springs) {
            int i = spring.first.first - offset;
            int j = spring.first.second - offset;
            if(i < 0 || i >= n_pro || j < 0 || j >= n_pro) {
                throw oxDNAException("The parameter file contains a spring between particles %d and %d, but only amino acids can be connected by springs", spring.first.first, spring.first.second);
            }
            per_residue[i].emplace_back(j, &spring.second);
            per_residue[j].emplace_back(i, &spring.second);
        }
        int n_entries = 2 * springs.size();

        _h_affected_indx = new int[n_pro + 1]();
        _h_affected = new int[n_entries]();
        _h_aff_gamma = new c_number[n_entries]();
        _h_aff_eqdist = new c_number[n_entries]();
        int idx = 0;
        for(int i = 0; i < n_pro; i++) {
            auto &neighs = per_residue[i];
            std::sort(neighs.begin(), neighs.end(), [](const std::pair<int, const ANMModel::Spring *> &a, const std::pair<int, const ANMModel::Spring *> &b) {
                return a.first < b.first;
            });
            _h_affected_indx[i] = idx;
            for(auto &neigh : neighs) {
                _h_affected[idx] = neigh.first;
                _h_aff_gamma[idx] = (c_number) neigh.second->k;
                _h_aff_eqdist[idx] = (c_number) neigh.second->r0;
                idx++;
            }
        }
        _h_affected_indx[n_pro] = idx;

        CUDA_SAFE_CALL(cudaMalloc(&_d_affected, n_entries * sizeof(int)));
        CUDA_SAFE_CALL(cudaMemcpy(_d_affected, _h_affected, n_entries * sizeof(int), cudaMemcpyHostToDevice));
        CUDA_SAFE_CALL(cudaMalloc(&_d_aff_gamma, n_entries * sizeof(c_number)));
        CUDA_SAFE_CALL(cudaMemcpy(_d_aff_gamma, _h_aff_gamma, n_entries * sizeof(c_number), cudaMemcpyHostToDevice));
        CUDA_SAFE_CALL(cudaMalloc(&_d_aff_eqdist, n_entries * sizeof(c_number)));
        CUDA_SAFE_CALL(cudaMemcpy(_d_aff_eqdist, _h_aff_eqdist, n_entries * sizeof(c_number), cudaMemcpyHostToDevice));
        CUDA_SAFE_CALL(cudaMalloc(&_d_affected_indx, (n_pro + 1) * sizeof(int)));
        CUDA_SAFE_CALL(cudaMemcpy(_d_affected_indx, _h_affected_indx, (n_pro + 1) * sizeof(int), cudaMemcpyHostToDevice));

        if(this->_anm.torsion()) {
            // bending/torsion parameters of the pair (i, i + 1), stored with the local index of i
            _ang_param_size = sizeof(c_number) * n_pro * 4;
            _h_ang_params = new c_number[n_pro * 4]();
            if(this->_anm.per_pair_stiffness()) {
                _h_ang_kbkt = new c_number[n_pro * 2]();
            }
            for(auto &ang : this->_anm.angular_parameters()) {
                int i = ang.first - offset;
                for(int k = 0; k < 4; k++) {
                    _h_ang_params[i * 4 + k] = (c_number) ang.second.cosines[k];
                }
                if(_h_ang_kbkt != nullptr) {
                    _h_ang_kbkt[i * 2] = (c_number) ang.second.k_bend;
                    _h_ang_kbkt[i * 2 + 1] = (c_number) ang.second.k_tor;
                }
            }
            CUDA_SAFE_CALL(cudaMalloc(&_d_ang_params, _ang_param_size));
            CUDA_SAFE_CALL(cudaMemcpy(_d_ang_params, _h_ang_params, _ang_param_size, cudaMemcpyHostToDevice));
            if(_h_ang_kbkt != nullptr) {
                CUDA_SAFE_CALL(cudaMalloc(&_d_ang_kbkt, sizeof(c_number) * n_pro * 2));
                CUDA_SAFE_CALL(cudaMemcpy(_d_ang_kbkt, _h_ang_kbkt, sizeof(c_number) * n_pro * 2, cudaMemcpyHostToDevice));
            }
        }

        OX_LOG(Logger::LOG_INFO, "Spring Parameters Size: %.2f MB", (n_entries * (sizeof(int) + 2 * sizeof(c_number)) + (n_pro + 1) * sizeof(int) + _ang_param_size) / (float) SQR(1024));
    }
    else OX_LOG(Logger::LOG_INFO, "Parfile: NONE, No protein parameters were filled");

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

    if(this->_anm.torsion()) {
        if (!this->_anm.per_pair_stiffness()) {
            _kbend = this->_anm.global_k_bend();
            _ktor = this->_anm.global_k_tor();

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
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(_ndna, &_n_nucleotides, sizeof(int)) );
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(_npro, &_n_amino_acids, sizeof(int)) );
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
