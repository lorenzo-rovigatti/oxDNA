/*
 * CUDALJInteraction.cu
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include "CUDALJInteraction.h"

#include "CUDA_LJ.cuh"
#include "../Lists/CUDASimpleVerletList.h"
#include "../Lists/CUDANoList.h"

CUDALJInteraction::CUDALJInteraction() {
	_edge_compatible = true;
}

CUDALJInteraction::~CUDALJInteraction() {

}

void CUDALJInteraction::get_settings(input_file &inp) {
	LJInteraction::get_settings(inp);
}

void CUDALJInteraction::cuda_init(int N) {
	CUDABaseInteraction::cuda_init(N);
	LJInteraction::init();

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_LJ_n, &this->_n, 3 * sizeof(int)));

	COPY_ARRAY_TO_CONSTANT(MD_sqr_rcut, this->_sqr_LJ_rcut, 3);
	COPY_ARRAY_TO_CONSTANT(MD_sqr_sigma, this->_sqr_sigma, 3);
	COPY_ARRAY_TO_CONSTANT(MD_epsilon, this->_epsilon, 3);
	COPY_ARRAY_TO_CONSTANT(MD_E_cut, this->_E_cut, 3);

	if(this->_use_edge) {
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_n_forces, &this->_n_forces, sizeof(int)));
	}
}

void CUDALJInteraction::compute_forces(CUDABaseList*lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox*d_box) {
	if(_use_edge) {
		lj_forces_edge
			<<<(lists->N_edges - 1)/(this->_launch_cfg.threads_per_block) + 1, this->_launch_cfg.threads_per_block>>>
			(d_poss, _d_edge_forces, lists->d_edge_list, lists->N_edges, d_box);

		this->_sum_edge_forces(d_forces);
		CUT_CHECK_ERROR("forces_second_step error -- lj");
	}
	else {
		lj_forces
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_forces, lists->d_matrix_neighs, lists->d_number_neighs, d_box);
		CUT_CHECK_ERROR("forces_second_step lj simple_lists error");
	}
}
