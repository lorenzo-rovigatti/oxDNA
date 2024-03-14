// Main interaction for CUDA Patchy Spring Particles


#include "CUDAPSP.cuh"

CUDA_PSPInteraction::CUDA_PSPInteraction():CUDABaseInteraction(), PSPInteraction() {

}

CUDA_PSPInteraction::~CUDA_PSPInteraction() {

}

void CUDA_PSPInteraction::get_settings(input_file &inp) {
    PSPInteraction::get_settings(inp);
}

void CUDA_PSPInteraction::cuda_init(int N) {
    CUDABaseInteraction::cuda_init(N);
    PSPInteraction::cuda_init(N);

    CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&d3bodyForce, N * sizeof(c_number4)));
    CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&d3bodyTorque, N * sizeof(c_number4)));

    
}

void CUDA_PSPInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torque, LR_bonds *d_bonds, CUDABox *d_box) {
    int N = CUDABaseInteraction::_N;
    if(_update_st) {
        CUDA_SAFE_CALL(cudaMemset(_d_st, 0, _N * sizeof(CUDAStressTensor)));
    }

    // PSP_forces
    //     <<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
    //     (d_poss, d_orientations, d_forces, d_torque, lists->d_matrix_neighs, lists->d_number_neighs, d_bonds, _rcut, _sqr_rcut, _update_st, _d_st, d_box);
}