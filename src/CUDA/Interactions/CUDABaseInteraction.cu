/*
 * CUDABaseInteraction.cu
 *
 *  Created on: 03/apr/2013
 *      Author: lorenzo
 */

#include "CUDABaseInteraction.h"

#include "../cuda_utils/CUDA_lr_common.cuh"
#include "../../Utilities/oxDNAException.h"

__global__ void sum_edge_forces_torques(c_number4 *edge_forces, c_number4 *forces, c_number4 *edge_torques, c_number4 *torques, int N, int n_forces) {
	if(IND >= N) return;

	c_number4 tot_force = forces[IND];
	c_number4 tot_torque = torques[IND];
	for(int i = 0; i < n_forces; i++) {
		c_number4 tmp_force = edge_forces[N * i + IND];
		tot_force.x += tmp_force.x;
		tot_force.y += tmp_force.y;
		tot_force.z += tmp_force.z;
		tot_force.w += tmp_force.w;
		edge_forces[N * i + IND] = make_c_number4(0, 0, 0, 0);

		c_number4 tmp_torque = edge_torques[N * i + IND];
		tot_torque.x += tmp_torque.x;
		tot_torque.y += tmp_torque.y;
		tot_torque.z += tmp_torque.z;
		tot_torque.w += tmp_torque.w;
		edge_torques[N * i + IND] = make_c_number4(0, 0, 0, 0);
	}

	forces[IND] = tot_force;
	torques[IND] = tot_torque;
}

__global__ void sum_edge_forces(c_number4 *edge_forces, c_number4 *forces, int N, int n_forces) {
	if(IND >= N) return;

	c_number4 tot_force = forces[IND];
	for(int i = 0; i < n_forces; i++) {
		c_number4 tmp_force = edge_forces[N * i + IND];
		tot_force.x += tmp_force.x;
		tot_force.y += tmp_force.y;
		tot_force.z += tmp_force.z;
		tot_force.w += tmp_force.w;
		edge_forces[N * i + IND] = make_c_number4(0, 0, 0, 0);
	}

	forces[IND] = tot_force;
}

CUDABaseInteraction::CUDABaseInteraction() {
	_use_edge = false;
	_n_forces = 1;

	_d_edge_forces = nullptr;
	_d_edge_torques = nullptr;
}

CUDABaseInteraction::~CUDABaseInteraction() {
	if(_use_edge) {
		if(_d_edge_forces != NULL) CUDA_SAFE_CALL(cudaFree(_d_edge_forces));
		if(_d_edge_torques != NULL) CUDA_SAFE_CALL(cudaFree(_d_edge_torques));
	}
}

void CUDABaseInteraction::_sum_edge_forces(c_number4 *d_forces) {
	sum_edge_forces
		<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
		(_d_edge_forces, d_forces, _N, _n_forces);
		CUT_CHECK_ERROR("sum_edge_forces error");
}

void CUDABaseInteraction::_sum_edge_forces_torques(c_number4 *d_forces, c_number4 *d_torques) {
	sum_edge_forces_torques
		<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
		(_d_edge_forces, d_forces, _d_edge_torques, d_torques, _N, _n_forces);
		CUT_CHECK_ERROR("sum_edge_forces_torques_error");
}

void CUDABaseInteraction::get_cuda_settings(input_file &inp) {
	int tmpi;
	if(getInputBoolAsInt(&inp, "use_edge", &tmpi, 0) == KEY_FOUND) {
		if(tmpi > 0) {
			_use_edge = true;
			getInputInt(&inp, "edge_n_forces", &_n_forces, 0);
			if(_n_forces < 1) throw oxDNAException("edge_n_forces must be > 0");
		}
	}
}

void CUDABaseInteraction::cuda_init(c_number box_side, int N) {
	_N = N;
	_box_side = box_side;

	if(_use_edge && _d_edge_forces == nullptr) {
		size_t size = sizeof(c_number4) * _N * _n_forces;
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<c_number4>(&_d_edge_forces, size));
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<c_number4>(&_d_edge_torques, size));

		CUDA_SAFE_CALL(cudaMemset(_d_edge_forces, 0, size));
		CUDA_SAFE_CALL(cudaMemset(_d_edge_torques, 0, size));
	}
}

void CUDABaseInteraction::set_launch_cfg(CUDA_kernel_cfg &launch_cfg) {
	_launch_cfg = launch_cfg;
}

void CUDABaseInteraction::_hb_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, float *hb_energies, int n_threads, bool *region_is_nearhb, CUDA_kernel_cfg _ffs_hb_precalc_kernel_cfg, CUDABox*d_box) {
	throw oxDNAException("On CUDA, FFS is only implemented for the DNA and RNA interactions");
}

void CUDABaseInteraction::_near_hb_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, bool *nearly_bonded_array, int n_threads, bool *region_is_nearhb, CUDA_kernel_cfg _ffs_hb_precalc_kernel_cfg, CUDABox*d_box) {
	throw oxDNAException("On CUDA, FFS is only implemented for the DNA and RNA interactions");
}

void CUDABaseInteraction::_dist_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, c_number *op_dists, int n_threads, CUDA_kernel_cfg _ffs_dist_precalc_kernel_cfg, CUDABox*d_box) {
	throw oxDNAException("On CUDA, FFS is only implemented for the DNA and RNA interactions");
}
