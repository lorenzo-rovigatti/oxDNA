/*
 * GpuUtils.cpp
 *
 *  Created on: 24/set/2010
 *      Author: lorenzo
 */

#include "CUDAUtils.h"
#include <curand_kernel.h>

size_t GpuUtils::_allocated_dev_mem = 0;

#ifndef OLD_ARCH
__global__ void print_array(int *v, int N) {
	for(int i = 0; i < N; i++) printf("%d %d\n", i, v[i]);
}

__global__ void print_array(float *v, int N) {
	for(int i = 0; i < N; i++) printf("%d %f\n", i, v[i]);
}

__global__ void print_array(double *v, int N) {
	for(int i = 0; i < N; i++) printf("%d %lf\n", i, v[i]);
}

__global__ void print_array(LR_double4 *v, int N) {
	for(int i = 0; i < N; i++) printf("%d %lf %lf %lf %lf\n", i, v[i].x, v[i].y, v[i].z, v[i].w);
}

__global__ void print_array(float4 *v, int N) {
	for(int i = 0; i < N; i++) printf("%d %lf %lf %lf %lf\n", i, v[i].x, v[i].y, v[i].z, v[i].w);
}

template<typename T>
__global__ void check_thresold(T *v, int N, int t) {
	for(int i = 0; i < N; i++) if(v[i] >= t) printf("%d %d\n", i, v[i]);
}

template<typename T>
void GpuUtils::print_device_array(T *v, int N) {
	print_array
		<<<1,1>>>
		(v, N);
	CUT_CHECK_ERROR("print_device_array error");
	cudaThreadSynchronize();
}

template<typename T>
void GpuUtils::check_device_thresold(T *v, int N, int t) {
	check_thresold<T>
		<<<1,1>>>
		(v, N, t);
	CUT_CHECK_ERROR("check_device_thresold error");
	cudaThreadSynchronize();
}
#endif

// template instantiation
template struct GPU_quat<float>;
template struct GPU_quat<double>;

template float GpuUtils::sum_4th_comp<float, float4>(float4 *, int);

#ifndef OLD_ARCH
template double GpuUtils::sum_4th_comp<double, double4>(double4 *, int);
template void GpuUtils::check_device_thresold<int>(int *, int, int);
template void GpuUtils::check_device_thresold<uint>(uint *, int, int);
template void GpuUtils::print_device_array<int>(int *v, int);
template void GpuUtils::print_device_array<float>(float *v, int);
template void GpuUtils::print_device_array<double>(double *v, int);
template void GpuUtils::print_device_array<LR_double4>(LR_double4 *v, int);
template void GpuUtils::print_device_array<float4>(float4 *v, int);
#endif
