/*
 * CUDABaseInteraction.cu
 *
 *  Created on: 03/apr/2013
 *      Author: lorenzo
 */

#include "CUDABaseInteraction.h"

#include "../cuda_utils/CUDA_lr_common.cuh"
#include "../../Utilities/oxDNAException.h"
#include "../../Utilities/ConfigInfo.h"

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/fill.h>
#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/transform_reduce.h>

__global__ void sum_edge_forces_torques(c_number4 __restrict__ *edge_forces, c_number4 __restrict__ *forces, c_number4 __restrict__ *edge_torques,
		c_number4 __restrict__ *torques, int N, int n_forces) {
	if(IND >= N) return;

	c_number4 tot_force = forces[IND];
	c_number4 tot_torque = torques[IND];
	for(int i = 0; i < n_forces; i++) {
		tot_force += edge_forces[N * i + IND];
		tot_torque += edge_torques[N * i + IND];
		edge_forces[N * i + IND] = make_c_number4(0.f, 0.f, 0.f, 0.f);
		edge_torques[N * i + IND] = make_c_number4(0.f, 0.f, 0.f, 0.f);
	}

	forces[IND] = tot_force;
	torques[IND] = tot_torque;
}

__global__ void sum_edge_forces(c_number4 __restrict__ *edge_forces, c_number4 __restrict__ *forces, int N, int n_forces) {
	if(IND >= N) return;

	c_number4 tot_force = forces[IND];
	for(int i = 0; i < n_forces; i++) {
		tot_force += edge_forces[N * i + IND];
		edge_forces[N * i + IND] = make_c_number4(0.f, 0.f, 0.f, 0.f);
	}

	forces[IND] = tot_force;
}

CUDABaseInteraction::CUDABaseInteraction() {

}

CUDABaseInteraction::~CUDABaseInteraction() {
	if(_use_edge) {
		if(_d_edge_forces != NULL) CUDA_SAFE_CALL(cudaFree(_d_edge_forces));
		if(_d_edge_torques != NULL) CUDA_SAFE_CALL(cudaFree(_d_edge_torques));
	}

	if(_d_st != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_st));
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
	int update_st_every = 0;
	getInputInt(&inp, "CUDA_update_stress_tensor_every", &update_st_every, 0);
	if(update_st_every > 0) {
		_update_st = true;
	}

	getInputBool(&inp, "use_edge", &_use_edge, 0);
	if(_use_edge) {
		if(!_edge_compatible) {
			throw oxDNAException("The selected CUDA interaction is not compatible with 'use_edge = true'");
		}

		getInputInt(&inp, "edge_n_forces", &_n_forces, 0);
		if(_n_forces < 1) {
			throw oxDNAException("edge_n_forces must be > 0");
		}
	}
}

void CUDABaseInteraction::cuda_init(int N) {
	_N = N;

	if(_use_edge && _d_edge_forces == nullptr) {
		size_t size = sizeof(c_number4) * _N * _n_forces;
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<c_number4>(&_d_edge_forces, size));
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<c_number4>(&_d_edge_torques, size));

		CUDA_SAFE_CALL(cudaMemset(_d_edge_forces, 0, size));
		CUDA_SAFE_CALL(cudaMemset(_d_edge_torques, 0, size));
	}

	if(_update_st) {
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_st, N * sizeof(CUDAStressTensor)));
	}
}

void CUDABaseInteraction::set_launch_cfg(CUDA_kernel_cfg &launch_cfg) {
	_launch_cfg = launch_cfg;
}

struct vel_to_st {
	__device__ CUDAStressTensor operator()(const c_number4 &v) const {
		return CUDAStressTensor(
			SQR(v.x),
			SQR(v.y),
			SQR(v.z),
			v.x * v.y,
			v.x * v.z,
			v.y * v.z
		);
	}
};

StressTensor CUDABaseInteraction::CPU_stress_tensor(c_number4 *vels) {
	thrust::device_ptr<CUDAStressTensor> t_st = thrust::device_pointer_cast(_d_st);
	CUDAStressTensor st_sum = thrust::reduce(t_st, t_st + _N, CUDAStressTensor());

	thrust::device_ptr<c_number4> t_vels = thrust::device_pointer_cast(vels);
	st_sum += thrust::transform_reduce(t_vels, t_vels + _N, vel_to_st(), CUDAStressTensor(), thrust::plus<CUDAStressTensor>());

	return st_sum.as_StressTensor();
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
