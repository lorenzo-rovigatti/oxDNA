/*
 * GpuUtils.cpp
 *
 *  Created on: 24/set/2010
 *      Author: lorenzo
 */

#include "CUDAUtils.h"
#include <curand_kernel.h>
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
#include <thrust/transform_reduce.h>

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

template<typename number4>
struct sum_number4 {
	__device__ number4 operator()(const number4& a, const number4& b) const {
		number4 res;
        	res.x = a.x + b.x;
		res.y = a.y + b.y;
		res.z = a.z + b.z;
        	res.w = a.w + b.w;
        	return res;
    	}
 };

template<typename number4>
struct number4_to_double {
	__device__ double operator()(const number4 &a) {
		return (double) a.w;
	}
};

template<typename number4>
number4 GpuUtils::sum_number4_on_GPU(number4 *dv, int N) {
	thrust::device_ptr<number4> t_dv = thrust::device_pointer_cast(dv);
	number4 zero = {0., 0., 0., 0.};
	return thrust::reduce(t_dv, t_dv + N, zero, sum_number4<number4>());
}

template<typename number4>
double GpuUtils::sum_number4_to_double_on_GPU(number4 *dv, int N) {
	thrust::device_ptr<number4> t_dv = thrust::device_pointer_cast(dv);
	return thrust::transform_reduce(t_dv, t_dv + N, number4_to_double<number4>(), 0., thrust::plus<double>());
}
// template instantiation
template struct GPU_quat<float>;
template struct GPU_quat<double>;

template float GpuUtils::sum_4th_comp<float, float4>(float4 *, int);
template float4 GpuUtils::sum_number4_on_GPU<float4>(float4 *dv, int N);
template LR_double4 GpuUtils::sum_number4_on_GPU<LR_double4>(LR_double4 *dv, int N);
template double GpuUtils::sum_number4_to_double_on_GPU<float4>(float4 *dv, int N);
template double GpuUtils::sum_number4_to_double_on_GPU<LR_double4>(LR_double4 *dv, int N);

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
