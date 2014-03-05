/*
 * CUDABaseThermostat.cu
 *
 *  Created on: Feb 15, 2013
 *      Author: rovigatti
 */

#include "CUDABaseThermostat.h"

template<typename number, typename number4>
CUDABaseThermostat<number, number4>::~CUDABaseThermostat() {
	CUDA_SAFE_CALL( cudaFree(_d_rand_state) );
}
template<typename number, typename number4>
void CUDABaseThermostat<number, number4>::get_cuda_settings(input_file &inp) {
	_launch_cfg.threads_per_block = 64;
	getInputInt(&inp, "threads_per_block", &_launch_cfg.threads_per_block, 0);
}

template<typename number, typename number4>
void CUDABaseThermostat<number, number4>::_setup_rand(int N) {
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<curandState>(&_d_rand_state, N * sizeof(curandState)) );

	_launch_cfg.blocks.x = N / _launch_cfg.threads_per_block + ((N % _launch_cfg.threads_per_block == 0) ? 0 : 1);
	if(_launch_cfg.blocks.x == 0) _launch_cfg.blocks.x = 1;
	_launch_cfg.blocks.y = _launch_cfg.blocks.z = 1;

	setup_curand<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>(_d_rand_state, _seed, N);
}

template class CUDABaseThermostat<float, float4>;
template class CUDABaseThermostat<double, LR_double4>;

