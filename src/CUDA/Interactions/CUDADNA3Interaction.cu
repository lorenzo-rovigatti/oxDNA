/*
 * CUDADNA3Interaction.cu
 *
 *  Created on: 13/may/25
 *      Author: lorenzo
 */

#include "CUDADNA3Interaction.h"

#include "CUDA_DNA3.cuh"
#include "../Lists/CUDASimpleVerletList.h"
#include "../Lists/CUDANoList.h"
#include "../../Interactions/DNA2Interaction.h"

#include "../CUDAUtils.h"

CUDADNA3Interaction::CUDADNA3Interaction() {
    _edge_compatible = true;
}

CUDADNA3Interaction::~CUDADNA3Interaction() {
    if(_d_is_strand_end != nullptr) {
        CUDA_SAFE_CALL(cudaFree(_d_is_strand_end));
    }
    if(_d_particle_types != nullptr) {
        CUDA_SAFE_CALL(cudaFree(_d_particle_types));
    }
}

void CUDADNA3Interaction::get_settings(input_file &inp) {
	Logger::instance()->disable_log("CUDADNA3Interaction");
    DNA3Interaction::get_settings(inp);
	Logger::instance()->enable_log("CUDADNA3Interaction");
}

std::vector<float> CUDADNA3Interaction::_convert_param_array(const MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> *src, int N_arrays) {
    std::vector<float> data(src[0].total_size * N_arrays);
    for(int i = 0; i < N_arrays; i++) {
        int offset = i * src[0].total_size;
        std::transform(src[i].data, src[i].data + src[i].total_size, data.begin() + offset, [](number v) { return static_cast<float>(v); });
    }
    return data;
}

void CUDADNA3Interaction::cuda_init(int N) {
    CUDABaseInteraction::cuda_init(N);
	Logger::instance()->disable_log("CUDADNA3Interaction");
    DNA3Interaction::init();
	Logger::instance()->enable_log("CUDADNA3Interaction");

    float f_copy = this->_hb_multiplier;
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_hb_multi, &f_copy, sizeof(float)));

    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));

    // oxDNA3

    // FENE
    std::vector<float> params = _convert_param_array(&this->_fene_r0_SD, 1);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_fene_r0_SD, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(&this->_fene_delta2_SD, 1);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_fene_delta2_SD, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(&this->_mbf_xmax_SD, 1);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_mbf_xmax_SD, params.data(), params.size() * sizeof(float)));
    COPY_NUMBER_TO_FLOAT(MD_fene_eps, this->_fene_eps);

    // excluded volume
    params = _convert_param_array(this->_excl_s, 7);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_excl_s, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(this->_excl_r, 7);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_excl_r, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(this->_excl_b, 7);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_excl_b, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(this->_excl_rc, 7);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_excl_rc, params.data(), params.size() * sizeof(float)));

    params = _convert_param_array(F1_SD_EPS, 2);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F1_SD_EPS, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F1_SD_A, 2);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F1_SD_A, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F1_SD_RC, 2);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F1_SD_RC, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F1_SD_R0, 2);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F1_SD_R0, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F1_SD_BLOW, 2);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F1_SD_BLOW, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F1_SD_BHIGH, 2);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F1_SD_BHIGH, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F1_SD_RLOW, 2);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F1_SD_RLOW, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F1_SD_RHIGH, 2);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F1_SD_RHIGH, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F1_SD_RCLOW, 2);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F1_SD_RCLOW, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F1_SD_RCHIGH, 2);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F1_SD_RCHIGH, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F1_SD_SHIFT, 2);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F1_SD_SHIFT, params.data(), params.size() * sizeof(float)));

    params = _convert_param_array(F2_SD_K, 4);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F2_SD_K, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F2_SD_RC, 4);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F2_SD_RC, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F2_SD_R0, 4);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F2_SD_R0, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F2_SD_BLOW, 4);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F2_SD_BLOW, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F2_SD_RLOW, 4);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F2_SD_RLOW, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F2_SD_RCLOW, 4);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F2_SD_RCLOW, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F2_SD_BHIGH, 4);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F2_SD_BHIGH, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F2_SD_RCHIGH, 4);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F2_SD_RCHIGH, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F2_SD_RHIGH, 4);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F2_SD_RHIGH, params.data(), params.size() * sizeof(float)));

    params = _convert_param_array(F4_SD_THETA_A, 21);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F4_SD_THETA_A, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F4_SD_THETA_B, 21);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F4_SD_THETA_B, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F4_SD_THETA_T0, 21);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F4_SD_THETA_T0, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F4_SD_THETA_TS, 21);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F4_SD_THETA_TS, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F4_SD_THETA_TC, 21);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F4_SD_THETA_TC, params.data(), params.size() * sizeof(float)));

    params = _convert_param_array(F5_SD_PHI_A, 4);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F5_SD_PHI_A, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F5_SD_PHI_B, 4);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F5_SD_PHI_B, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F5_SD_PHI_XC, 4);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F5_SD_PHI_XC, params.data(), params.size() * sizeof(float)));
    params = _convert_param_array(F5_SD_PHI_XS, 4);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_F5_SD_PHI_XS, params.data(), params.size() * sizeof(float)));

    // end oxDNA3

    std::vector<uint8_t> particle_types(this->_N);
    for(int i = 0; i < _N; i++) {
        particle_types[i] = CONFIG_INFO->particles()[i]->type;
    }
    CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<uint8_t>(&_d_particle_types, sizeof(uint8_t) * _N));
    CUDA_SAFE_CALL(cudaMemcpy(_d_particle_types, particle_types.data(), sizeof(uint8_t) * _N, cudaMemcpyHostToDevice));

    if(this->_use_edge) CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_n_forces, &this->_n_forces, sizeof(int)));

    COPY_NUMBER_TO_FLOAT(MD_dh_RC, _debye_huckel_RC);
    COPY_NUMBER_TO_FLOAT(MD_dh_RHIGH, _debye_huckel_RHIGH);
    COPY_NUMBER_TO_FLOAT(MD_dh_prefactor, _debye_huckel_prefactor);
    COPY_NUMBER_TO_FLOAT(MD_dh_B, _debye_huckel_B);
    COPY_NUMBER_TO_FLOAT(MD_dh_minus_kappa, _minus_kappa);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_dh_half_charged_ends, &_debye_huckel_half_charged_ends, sizeof(bool)));
}

void CUDADNA3Interaction::_on_T_update() {
    cuda_init(_N);
}

void CUDADNA3Interaction::_init_strand_ends(LR_bonds *d_bonds) {
    CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&_d_is_strand_end, sizeof(int) * _N));
    init_DNA3_strand_ends<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>(_d_is_strand_end, d_bonds, _N);
}

void CUDADNA3Interaction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox*d_box) {
    if(_d_is_strand_end == nullptr) {
        _init_strand_ends(d_bonds);
    }

    if(_update_st) {
        CUDA_SAFE_CALL(cudaMemset(_d_st, 0, _N * sizeof(CUDAStressTensor)));
    }

    if(_use_edge) {
        if(_n_forces == 1) { // we can directly use d_forces and d_torques so that no sum is required
            DNA3_forces_edge_nonbonded
                <<<(lists->N_edges - 1)/(_launch_cfg.threads_per_block) + 1, _launch_cfg.threads_per_block>>>
                (d_poss, d_orientations, d_forces, d_torques, lists->d_edge_list, lists->N_edges, d_bonds,
                _d_particle_types, _update_st, _d_st, d_box);
        }
        else { // sum required, somewhat slower
            DNA3_forces_edge_nonbonded
                <<<(lists->N_edges - 1)/(_launch_cfg.threads_per_block) + 1, _launch_cfg.threads_per_block>>>
                (d_poss, d_orientations, _d_edge_forces, _d_edge_torques, lists->d_edge_list, lists->N_edges, 
                d_bonds, _d_particle_types, _update_st, _d_st, d_box);

            _sum_edge_forces_torques(d_forces, d_torques);
        }

        DNA3_forces_edge_bonded
            <<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
            (d_poss, d_orientations, d_forces, d_torques, d_bonds, _use_mbf, _mbf_finf, _d_particle_types, _update_st, _d_st);
    }
    else {
        DNA3_forces
            <<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
            (d_poss, d_orientations, d_forces, d_torques, lists->d_matrix_neighs, lists->d_number_neighs,
            d_bonds, _use_mbf, _mbf_finf, _d_particle_types, _update_st, _d_st, d_box);
        CUT_CHECK_ERROR("forces_second_step simple_lists error");
    }

    // GpuUtils::print_device_array(d_forces, this->_N);
    // GpuUtils::print_device_array(d_torques, this->_N);
}
