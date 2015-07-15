/*
 * MD_CUDAMixedBackend.cu
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo */

#include "MD_CUDAMixedBackend.h"

#include "CUDA_mixed.cuh"

CUDAMixedBackend::CUDAMixedBackend() : MD_CUDABackend<float, float4>() {
	_d_possd = NULL;
	_d_velsd = NULL;
	_d_Lsd = NULL;
	_d_orientationsd = NULL;
}

CUDAMixedBackend::~CUDAMixedBackend(){
	CUDA_SAFE_CALL( cudaFree(_d_possd) );
	CUDA_SAFE_CALL( cudaFree(_d_orientationsd) );
	CUDA_SAFE_CALL( cudaFree(_d_velsd) );
	CUDA_SAFE_CALL( cudaFree(_d_Lsd) );
}

void CUDAMixedBackend::init() {
	MD_CUDABackend<float, float4>::init();

	_vec_sized = ((size_t)_N) * sizeof(LR_double4);
	_orient_sized = ((size_t)_N) * sizeof(GPU_quat<double>);
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<LR_double4>(&_d_possd, _vec_sized) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<LR_double4>(&_d_velsd, _vec_sized) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<LR_double4>(&_d_Lsd, _vec_sized) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<GPU_quat<double> >(&_d_orientationsd, _orient_sized) );

	_float4_to_LR_double4(_d_poss, _d_possd);
	_quat_float_to_quat_double(_d_orientations, _d_orientationsd); 
	_float4_to_LR_double4(_d_vels, _d_velsd);
	_float4_to_LR_double4(_d_Ls, _d_Lsd);
}

void CUDAMixedBackend::_init_CUDA_MD_symbols() {
	MD_CUDABackend<float, float4>::_init_CUDA_MD_symbols();
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_sqr_verlet_skin, &_sqr_verlet_skin, sizeof(float)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dt, &_dt, sizeof(float)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &_N, sizeof(int)) );
}

void CUDAMixedBackend::_float4_to_LR_double4(float4 *src, LR_double4 *dest) {
	float4_to_LR_double4
		<<<_particles_kernel_cfg.blocks, _particles_kernel_cfg.threads_per_block>>>
		(src, dest);
	CUT_CHECK_ERROR("float4_to_LR_double4 error");
}

void CUDAMixedBackend::_LR_double4_to_float4(LR_double4 *src, float4 *dest) {
	LR_double4_to_float4
		<<<_particles_kernel_cfg.blocks, _particles_kernel_cfg.threads_per_block>>>
		(src, dest);
	CUT_CHECK_ERROR("LR_double4_to_float4 error");
}

void CUDAMixedBackend::_quat_float_to_quat_double(GPU_quat<float> *src, GPU_quat<double> *dest) {
	float4_to_LR_double4
		<<<_particles_kernel_cfg.blocks, _particles_kernel_cfg.threads_per_block>>>
		(src, dest);
	CUT_CHECK_ERROR("LR_matrix_float_to_matrix_double error");
}

void CUDAMixedBackend::_quat_double_to_quat_float(GPU_quat<double> *src, GPU_quat<float> *dest) {
	LR_double4_to_float4
		<<<_particles_kernel_cfg.blocks, _particles_kernel_cfg.threads_per_block>>>
		(src, dest);
	CUT_CHECK_ERROR("LR_matrix_double_to_matrix_float error");
}

void CUDAMixedBackend::_first_step() {
	first_step_mixed
		<<<_particles_kernel_cfg.blocks, _particles_kernel_cfg.threads_per_block>>>
		(_d_poss, _d_orientations, _d_possd, _d_orientationsd, _d_list_poss, _d_velsd, _d_Lsd, _d_forces, _d_torques, _d_are_lists_old, this->_any_rigid_body);
}

void CUDAMixedBackend::_sort_particles() {
	_LR_double4_to_float4(_d_possd, _d_poss);
	_LR_double4_to_float4(_d_velsd, _d_vels);
	_LR_double4_to_float4(_d_Lsd, _d_Ls);
	_quat_double_to_quat_float(_d_orientationsd, _d_orientations);
	MD_CUDABackend<float, float4>::_sort_particles();
	_quat_float_to_quat_double(_d_orientations, _d_orientationsd);
	_float4_to_LR_double4(_d_poss, _d_possd);
	_float4_to_LR_double4(_d_vels, _d_velsd);
	_float4_to_LR_double4(_d_Ls, _d_Lsd);
}

void CUDAMixedBackend::_forces_second_step() {
	this->_set_external_forces();
	_cuda_interaction->compute_forces(_cuda_lists, _d_poss, _d_orientations, _d_forces, _d_torques, _d_bonds);

	second_step_mixed
		<<<_particles_kernel_cfg.blocks, _particles_kernel_cfg.threads_per_block>>>
		(_d_velsd, _d_Lsd, _d_forces, _d_torques, this->_any_rigid_body);
	CUT_CHECK_ERROR("second_step_mixed");
}

void CUDAMixedBackend::_gpu_to_host_particles() {
	// probably useless
	if(_d_possd != NULL) {
		_LR_double4_to_float4(_d_possd, _d_poss);
		_quat_double_to_quat_float(_d_orientationsd, _d_orientations);
		_LR_double4_to_float4(_d_velsd, _d_vels);
		_LR_double4_to_float4(_d_Lsd, _d_Ls);
	}

	MD_CUDABackend<float, float4>::_gpu_to_host_particles();
}

void CUDAMixedBackend::_host_particles_to_gpu() {
	MD_CUDABackend<float, float4>::_host_particles_to_gpu();

	// the first time this method gets called all these arrays have not been
	// allocated yet. It's a bit of a hack but it's needed
	if(_d_possd != NULL) {
		_float4_to_LR_double4(_d_poss, _d_possd);
		_quat_float_to_quat_double(_d_orientations, _d_orientationsd);
		_float4_to_LR_double4(_d_vels, _d_velsd);
		_float4_to_LR_double4(_d_Ls, _d_Lsd);
	}
}

void CUDAMixedBackend::_thermalize(llint curr_step) {
	if(_cuda_thermostat->would_activate(curr_step)) {
		_LR_double4_to_float4(_d_velsd, _d_vels);
		_LR_double4_to_float4(_d_Lsd, _d_Ls);
		MD_CUDABackend<float, float4>::_thermalize(curr_step);
		_float4_to_LR_double4(_d_vels, _d_velsd);
		_float4_to_LR_double4(_d_Ls, _d_Lsd);
	}
}
