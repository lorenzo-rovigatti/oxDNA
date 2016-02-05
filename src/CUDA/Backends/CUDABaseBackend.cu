/*
 * CUDABaseBackend.cpp
 *
 *  Created on: 25/nov/2010
 *      Author: lorenzo
 */

#include <thrust/sort.h>
#include <thrust/device_ptr.h>

#include "CUDABaseBackend.h"
#include "../Lists/CUDAListFactory.h"
#include "../Interactions/CUDAInteractionFactory.h"
#include "../../Utilities/oxDNAException.h"

using namespace std;

template<typename number, typename number4>
CUDABaseBackend<number, number4>::CUDABaseBackend() : _device_number(0), _sort_every(0) {
	_particles_kernel_cfg.blocks = dim3(1, 1, 1);
	_particles_kernel_cfg.threads_per_block = 0;
	_particles_kernel_cfg.shared_mem = 0;

	_device_number = -1;
	_sqr_verlet_skin = 0.f;

	_cuda_lists = NULL;
	_cuda_interaction = NULL;
	_d_poss = NULL;
	_d_bonds = NULL;
	_d_orientations = NULL;
	_d_list_poss = NULL;
	_d_are_lists_old = NULL;
	_d_hindex = NULL;
	_d_sorted_hindex = NULL;
	_d_inv_sorted_hindex = NULL;
	_d_buff_poss = NULL;
	_d_buff_bonds = NULL;
	_d_buff_orientations = NULL;
	_h_poss = NULL;
	_h_orientations = NULL;
	_h_bonds = NULL;
}

template<typename number, typename number4>
CUDABaseBackend<number, number4>::~CUDABaseBackend() {
	if (_cuda_lists != NULL) {
		_cuda_lists->clean();
		delete _cuda_lists;
	}
	if (_cuda_interaction != NULL) delete _cuda_interaction;

	if (_d_poss != NULL){
		CUDA_SAFE_CALL( cudaFree(_d_poss) );
		CUDA_SAFE_CALL( cudaFree(_d_bonds) );
		CUDA_SAFE_CALL( cudaFree(_d_orientations) );
		CUDA_SAFE_CALL( cudaFree(_d_list_poss) );
		CUDA_SAFE_CALL( cudaFreeHost(_d_are_lists_old) );
	}

	if(_sort_every > 0) {
		if (_d_hindex != NULL){
			CUDA_SAFE_CALL( cudaFree(_d_hindex) );
			CUDA_SAFE_CALL( cudaFree(_d_sorted_hindex) );
			CUDA_SAFE_CALL( cudaFree(_d_inv_sorted_hindex) );
			CUDA_SAFE_CALL( cudaFree(_d_buff_poss) );
			CUDA_SAFE_CALL( cudaFree(_d_buff_bonds) );
			CUDA_SAFE_CALL( cudaFree(_d_buff_orientations) );
		}
	}

	if (_h_poss != NULL) delete[] _h_poss;
	if (_h_orientations != NULL) delete[] _h_orientations;
	if (_h_bonds != NULL) delete[] _h_bonds;
}

template<typename number, typename number4>
void CUDABaseBackend<number, number4>::_host_to_gpu() {
	CUDA_SAFE_CALL( cudaMemcpy(_d_poss, _h_poss, _vec_size, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_bonds, _h_bonds, _bonds_size, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_orientations, _h_orientations, _orient_size, cudaMemcpyHostToDevice) );
}

template<typename number, typename number4>
void CUDABaseBackend<number, number4>::_gpu_to_host() {
	CUDA_SAFE_CALL( cudaMemcpy(_h_poss, _d_poss, _vec_size, cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(_h_bonds, _d_bonds, _bonds_size, cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(_h_orientations, _d_orientations, _orient_size, cudaMemcpyDeviceToHost) );
}

template<typename number, typename number4>
void CUDABaseBackend<number, number4>::get_settings(input_file &inp) {
	if(getInputInt(&inp, "CUDA_device", &_device_number, 0) == KEY_NOT_FOUND) {
		OX_LOG(Logger::LOG_INFO, "CUDA device not specified");
		_device_number = -1;
	}
	else OX_LOG(Logger::LOG_INFO, "Using CUDA device %d", _device_number);

	if(getInputInt(&inp, "CUDA_sort_every", &_sort_every, 0) == KEY_NOT_FOUND)
		OX_LOG(Logger::LOG_INFO, "CUDA sort_every not specified, using 0");

	getInputInt(&inp, "threads_per_block", &_particles_kernel_cfg.threads_per_block, 0);

	float verlet_skin;
	if(getInputFloat(&inp, "verlet_skin", &verlet_skin, 0) == KEY_FOUND) _sqr_verlet_skin = SQR(verlet_skin);

	_cuda_interaction = CUDAInteractionFactory::make_interaction<number, number4>(inp);
	_cuda_interaction->get_settings(inp);
	_cuda_interaction->get_cuda_settings(inp);

	_cuda_lists = CUDAListFactory::make_list<number, number4>(inp);
	_cuda_lists->get_settings(inp);

	// check that the box is cubic
	string my_box;
	if(getInputString(&inp, "box_type", my_box, 0) == KEY_FOUND) if(my_box != "cubic") throw oxDNAException("The CUDA backend only supports cubic boxes");

	string reload_from;
	if(getInputString(&inp, "reload_from", reload_from, 0) == KEY_FOUND) throw oxDNAException("The CUDA backend does not support reloading checkpoints, owing to its intrisincally stochastic nature");
}

template<typename number, typename number4>
void CUDABaseBackend<number, number4>::_choose_device () {
	OX_LOG(Logger::LOG_INFO, "Choosing device automatically");

	int ndev = -1, trydev = 0;
	cudaDeviceProp tryprop;

	cudaGetDeviceCount (&ndev);
	OX_LOG(Logger::LOG_INFO, "Computer has %i devices", ndev);
	while (trydev < ndev) {
		OX_LOG(Logger::LOG_INFO, " - Trying device %i", trydev);
		tryprop = get_device_prop (trydev);
		OX_LOG(Logger::LOG_INFO, " -- device %i has properties %i.%i", trydev, tryprop.major, tryprop.minor);
		if (tryprop.major < 2 && tryprop.minor <= 2) {
			OX_LOG(Logger::LOG_INFO, " -- Device properties are not good. Skipping it", trydev);
			trydev ++;
			continue;
		}
		set_device (trydev);
		int *dummyptr = NULL;
		cudaError_t ggg = GpuUtils::LR_cudaMalloc<int> (&dummyptr, (size_t)sizeof(int));
		if(ggg == cudaSuccess) {
			OX_LOG(Logger::LOG_INFO, " -- using device %i", trydev);
			cudaFree (dummyptr);
			break;
		}
		else {
			OX_LOG(Logger::LOG_INFO, " -- device %i not available ...", trydev);
		}
		trydev ++;
	}

	if (trydev == ndev) throw oxDNAException("No suitable devices available");

	OX_LOG(Logger::LOG_INFO, " --- Running on device %i", trydev);
	_device_prop = get_device_prop(trydev);
	_device_number = trydev;
	// gpu device chosen
}

template<typename number, typename number4>
void CUDABaseBackend<number, number4>::init_cuda(ConfigInfo<number> &config_info) {
	if(_device_number < 0)	_choose_device();
	set_device(_device_number);
	_device_prop = get_device_prop(_device_number);

	CUDA_SAFE_CALL( cudaThreadSetCacheConfig(cudaFuncCachePreferL1) );

	_prv_config_info = config_info;
	number box_side = config_info.box->box_sides().x;
	int N = *config_info.N;

	_cuda_interaction->cuda_init(box_side, N);

	_vec_size = sizeof(number4) * N;
	_orient_size = sizeof(GPU_quat<number>) * N;
	_bonds_size = sizeof(LR_bonds) * N;

	// GPU memory allocations
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number4>(&_d_poss, _vec_size) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<LR_bonds>(&_d_bonds, _bonds_size) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<GPU_quat<number>  >(&_d_orientations, _orient_size) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number4>(&_d_list_poss, _vec_size) );
	CUDA_SAFE_CALL( cudaMallocHost(&_d_are_lists_old, sizeof(bool), cudaHostAllocDefault) );

	CUDA_SAFE_CALL( cudaMemset(_d_list_poss, 0, _vec_size) );

	// CPU memory allocations
	_h_poss = new number4[N];
	_h_orientations = new GPU_quat<number>[N];
	_h_bonds = new LR_bonds[N];

	// setup kernels' configurations
	_init_CUDA_kernel_cfgs();
	_cuda_lists->init(N, box_side, _cuda_interaction->get_cuda_rcut());

	if(_sort_every > 0) {
		int uns = 0;

		// fixed value for depth (8): changing this value does not significantly affect performances
		init_hilb_symbols(N, uns, 8, (float) box_side);

		CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&_d_hindex, N*sizeof(int)) );
		CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&_d_sorted_hindex, N*sizeof(int)) );
		CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&_d_inv_sorted_hindex, N*sizeof(int)) );
		CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number4>(&_d_buff_poss, _vec_size) );
		CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<LR_bonds>(&_d_buff_bonds, _bonds_size) );
		CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<GPU_quat<number>  >(&_d_buff_orientations, _orient_size) );

		reset_sorted_hindex
			<<<_particles_kernel_cfg.blocks, _particles_kernel_cfg.threads_per_block>>>
			(_d_sorted_hindex);
	}
}

template<typename number, typename number4>
void CUDABaseBackend<number, number4>::_init_CUDA_kernel_cfgs() {
	if(_particles_kernel_cfg.threads_per_block == 0) {
		_particles_kernel_cfg.threads_per_block = 2*_device_prop.warpSize;
		OX_LOG(Logger::LOG_INFO, "threads_per_block was not specified or set to 0. The default value (%d) will be used", 2*_device_prop.warpSize);
	}

	int N = *_prv_config_info.N;
	_particles_kernel_cfg.blocks.x = N / _particles_kernel_cfg.threads_per_block + ((N % _particles_kernel_cfg.threads_per_block == 0) ? 0 : 1);
	if(_particles_kernel_cfg.blocks.x == 0) _particles_kernel_cfg.blocks.x = 1;
	_particles_kernel_cfg.blocks.y = _particles_kernel_cfg.blocks.z = 1;

	_cuda_interaction->set_launch_cfg(_particles_kernel_cfg);

	OX_DEBUG("Particle kernel cfg: threads_per_block = %d, blocks = (%d, %d, %d)", _particles_kernel_cfg.threads_per_block,
			_particles_kernel_cfg.blocks.x, _particles_kernel_cfg.blocks.y, _particles_kernel_cfg.blocks.z);
}

template<typename number, typename number4>
void CUDABaseBackend<number, number4>::_sort_index() {
	reset_sorted_hindex
		<<<_particles_kernel_cfg.blocks, _particles_kernel_cfg.threads_per_block>>>
		(_d_sorted_hindex);
	CUT_CHECK_ERROR("reset_sorted_hindex error");

	hilbert_curve<number4>
		<<<_particles_kernel_cfg.blocks, _particles_kernel_cfg.threads_per_block>>>
		(_d_poss, _d_hindex);
	CUT_CHECK_ERROR("hilbert_curve error");

	thrust::device_ptr<int> _d_hindex_p(_d_hindex);
	thrust::device_ptr<int> _d_sorted_hindex_p(_d_sorted_hindex);
	// sort d_sorted_hindex by using d_hindex
	thrust::sort_by_key(_d_hindex_p, _d_hindex_p + *_prv_config_info.N, _d_sorted_hindex_p);
	get_inverted_sorted_hindex
		<<<_particles_kernel_cfg.blocks, _particles_kernel_cfg.threads_per_block>>>
		(_d_sorted_hindex, _d_inv_sorted_hindex);
}

// template instantiations
template class CUDABaseBackend<float, float4>;
template class CUDABaseBackend<double, LR_double4>;
