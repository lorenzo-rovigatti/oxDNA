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

template<typename number, typename number4>
CUDALJInteraction<number, number4>::CUDALJInteraction() {

}

template<typename number, typename number4>
CUDALJInteraction<number, number4>::~CUDALJInteraction() {

}

template<typename number, typename number4>
void CUDALJInteraction<number, number4>::get_settings(input_file &inp) {
	LJInteraction<number>::get_settings(inp);
}

template<typename number, typename number4>
void CUDALJInteraction<number, number4>::cuda_init(number box_side, int N) {
	CUDABaseInteraction<number, number4>::cuda_init(box_side, N);
	LJInteraction<number>::init();

	float f_copy = box_side;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_box_side, &f_copy, sizeof(float)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &N, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_LJ_n, &this->_n, 3*sizeof(int)) );

	COPY_ARRAY_TO_CONSTANT(MD_sqr_rcut, this->_sqr_LJ_rcut, 3);
	COPY_ARRAY_TO_CONSTANT(MD_sqr_sigma, this->_sqr_sigma, 3);
	COPY_ARRAY_TO_CONSTANT(MD_epsilon, this->_epsilon, 3);

	if(this->_use_edge) CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_n_forces, &this->_n_forces, sizeof(int)) );
}

template<typename number, typename number4>
void CUDALJInteraction<number, number4>::compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds) {
	CUDASimpleVerletList<number, number4> *_v_lists = dynamic_cast<CUDASimpleVerletList<number, number4> *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) {
				lj_forces_edge<number, number4>
					<<<(_v_lists->_N_edges - 1)/(this->_launch_cfg.threads_per_block) + 1, this->_launch_cfg.threads_per_block>>>
					//(d_poss, d_forces, d_torques, _v_lists->_d_edge_list, _v_lists->_N_edges);
					(d_poss, this->_d_edge_forces, _v_lists->_d_edge_list, _v_lists->_N_edges);

				this->_sum_edge_forces(d_forces);

				// potential for removal here
				cudaThreadSynchronize();
				CUT_CHECK_ERROR("forces_second_step error -- lj");
			}
			else {
				lj_forces<number, number4>
					<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
					(d_poss, d_forces, _v_lists->_d_matrix_neighs, _v_lists->_d_number_neighs);
				CUT_CHECK_ERROR("forces_second_step lj simple_lists error");
			}
	}

	CUDANoList<number, number4> *_no_lists = dynamic_cast<CUDANoList<number, number4> *>(lists);
	if(_no_lists != NULL) {
		lj_forces<number, number4>
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_forces);
		CUT_CHECK_ERROR("forces_second_step lj no_lists error");
	}
}

template class CUDALJInteraction<float, float4>;
template class CUDALJInteraction<double, LR_double4>;
