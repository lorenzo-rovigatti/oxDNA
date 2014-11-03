/*
 * CUDABaseInteraction.cu
 *
 *  Created on: 03/apr/2013
 *      Author: lorenzo
 */

#include "CUDABaseInteraction.h"

#include "../cuda_utils/CUDA_lr_common.cuh"
#include "../../Utilities/oxDNAException.h"

template <typename number, typename number4>
__global__ void sum_edge_forces_torques(number4 *edge_forces, number4 *forces, number4 *edge_torques, number4 *torques, int N, int n_forces) {
	if(IND >= N) return;

	number4 tot_force = forces[IND];
	number4 tot_torque = torques[IND];
	for(int i = 0; i < n_forces; i++) {
		number4 tmp_force = edge_forces[N*i + IND];
		tot_force.x += tmp_force.x;
		tot_force.y += tmp_force.y;
		tot_force.z += tmp_force.z;
		tot_force.w += tmp_force.w;
		edge_forces[N*i + IND] = make_number4<number, number4>(0, 0, 0, 0);

		number4 tmp_torque = edge_torques[N*i + IND];
		tot_torque.x += tmp_torque.x;
		tot_torque.y += tmp_torque.y;
		tot_torque.z += tmp_torque.z;
		tot_torque.w += tmp_torque.w;
		edge_torques[N*i + IND] = make_number4<number, number4>(0, 0, 0, 0);
	}

	forces[IND] = tot_force;
	torques[IND] = tot_torque;
}

template <typename number, typename number4>
__global__ void sum_edge_forces(number4 *edge_forces, number4 *forces, int N, int n_forces) {
	if(IND >= N) return;

	number4 tot_force = forces[IND];
	for(int i = 0; i < n_forces; i++) {
		number4 tmp_force = edge_forces[N*i + IND];
		tot_force.x += tmp_force.x;
		tot_force.y += tmp_force.y;
		tot_force.z += tmp_force.z;
		tot_force.w += tmp_force.w;
		edge_forces[N*i + IND] = make_number4<number, number4>(0, 0, 0, 0);
	}

	forces[IND] = tot_force;
}

template<typename number, typename number4>
CUDABaseInteraction<number, number4>::CUDABaseInteraction() {
	_use_edge = false;
	_n_forces = 1;

	_d_edge_forces = NULL;
	_d_edge_torques = NULL;
}

template<typename number, typename number4>
CUDABaseInteraction<number, number4>::~CUDABaseInteraction() {
	if(_use_edge) {
		if (_d_edge_forces != NULL) CUDA_SAFE_CALL( cudaFree(_d_edge_forces) );
		if (_d_edge_torques != NULL) CUDA_SAFE_CALL( cudaFree(_d_edge_torques) );
	}
}

template<typename number, typename number4>
void CUDABaseInteraction<number, number4>::_sum_edge_forces(number4 *d_forces) {
	sum_edge_forces<number, number4>
		<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
		(_d_edge_forces, d_forces, _N, _n_forces);
	CUT_CHECK_ERROR("sum_edge_forces error");
}

template<typename number, typename number4>
void CUDABaseInteraction<number, number4>::_sum_edge_forces_torques(number4 *d_forces, number4 *d_torques) {
	sum_edge_forces_torques<number, number4>
		<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
		(_d_edge_forces, d_forces, _d_edge_torques, d_torques, _N, _n_forces);
	CUT_CHECK_ERROR("sum_edge_forces_torques_error");
}

template<typename number, typename number4>
void CUDABaseInteraction<number, number4>::get_cuda_settings(input_file &inp) {
	int tmpi;
	if(getInputBoolAsInt(&inp, "use_edge", &tmpi, 0) == KEY_FOUND) {
		if(tmpi > 0) {
			_use_edge = true;
			getInputInt(&inp, "edge_n_forces", &_n_forces, 0);
			if(_n_forces < 1) throw oxDNAException("edge_n_forces must be > 0");
		}
	}
}

template<typename number, typename number4>
void CUDABaseInteraction<number, number4>::cuda_init(number box_side, int N) {
	_N = N;
	_box_side = box_side;

	if(_use_edge) {
		size_t size = sizeof(number4)*_N*_n_forces;
		CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number4>(&_d_edge_forces, size) );
		CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number4>(&_d_edge_torques, size) );

		CUDA_SAFE_CALL( cudaMemset(_d_edge_forces, 0, size) );
		CUDA_SAFE_CALL( cudaMemset(_d_edge_torques, 0, size) );
	}
}

template<typename number, typename number4>
void CUDABaseInteraction<number, number4>::set_launch_cfg(CUDA_kernel_cfg &launch_cfg) {
	_launch_cfg = launch_cfg;
}

template<typename number, typename number4>
void CUDABaseInteraction<number, number4>::_hb_op_precalc(number4 *poss, GPU_quat<number> *orientations, int *op_pairs1, int *op_pairs2, float *hb_energies, int n_threads, bool *region_is_nearhb, CUDA_kernel_cfg _ffs_hb_precalc_kernel_cfg) {
	throw oxDNAException ("On CUDA, FFS is only implemented for the DNA and RNA interactions");
}

template<typename number, typename number4>
void CUDABaseInteraction<number, number4>::_near_hb_op_precalc(number4 *poss, GPU_quat<number> *orientations, int *op_pairs1, int *op_pairs2, bool *nearly_bonded_array, int n_threads, bool *region_is_nearhb, CUDA_kernel_cfg  _ffs_hb_precalc_kernel_cfg) {
	throw oxDNAException ("On CUDA, FFS is only implemented for the DNA and RNA interactions");
}

template<typename number, typename number4>
void CUDABaseInteraction<number, number4>::_dist_op_precalc(number4 *poss, GPU_quat<number> *orientations, int *op_pairs1, int *op_pairs2, number *op_dists, int n_threads, CUDA_kernel_cfg _ffs_dist_precalc_kernel_cfg) {
	throw oxDNAException ("On CUDA, FFS is only implemented for the DNA and RNA interactions");
}

template class CUDABaseInteraction<float, float4>;
template class CUDABaseInteraction<double, LR_double4>;
