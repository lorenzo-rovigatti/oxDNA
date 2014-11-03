/**
 * @file    CUDAUtils.h
 * @date    24/set/2010
 * @author  lorenzo
 *
 *
 */

#ifndef GPUUTILS_H_
#define GPUUTILS_H_

#include <cuda.h>
#include <cuda_runtime_api.h>
#include "../defs.h"
#include "cuda_utils/helper_cuda.h"

/// CUDA_SAFE_CALL replacement for backwards compatibility (CUDA < 5.0)
#define CUDA_SAFE_CALL(x) checkCudaErrors(x);
/// CUT_CHECK_ERROR replacement for backwards compatibility (CUDA < 5.0)
#define CUT_CHECK_ERROR(x) getLastCudaError(x);

/// threads per block
#define TINBLOCK (blockDim.x*blockDim.y)
/// number of blocks
#define NBLOCKS (gridDim.x*gridDim.y)
/// number of threads
#define NTHREADS (NBLOCKS * TINBLOCK)

/// thread id relative to its block
#define TID (blockDim.x*threadIdx.y + threadIdx.x)
/// block id
#define BID (gridDim.x*blockIdx.y + blockIdx.x)
/// thread id
#define IND (TINBLOCK * BID + TID)

#define CUDA_LRACOS(x) (((x) >= (number)1) ? (number) 0 : ((x) <= (number)-1) ? (number) PI : acosf(x))
#define CUDA_DOT(a, b) (a.x*b.x + a.y*b.y + a.z*b.z)

#define COPY_ARRAY_TO_CONSTANT(dest, src, size) {\
		float *val = new float[(size)];\
		for(int i = 0; i < (size); i++) val[i] = (float) ((src)[i]);\
		CUDA_SAFE_CALL(cudaMemcpyToSymbol((dest), val, (size)*sizeof(float)))\
		delete[] val; }

/**
 * @brief Utility struct used by CUDA class to store information about kernel configurations.
 */
typedef struct CUDA_kernel_cfg {
	dim3 blocks;
	int threads_per_block;
	int shared_mem;
} CUDA_kernel_cfg;

/**
 * @brief We need this struct because the fourth element of such a structure must be a float or _float_as_int will not work.
 */
typedef struct __align__(16) {
	double x, y, z;
	float w;
} LR_double4;

/**
 * @brief It keeps track of neighbours along 3" and 5" directions.
 */
typedef struct __align__(8) {
	int n3, n5;
} LR_bonds;

/*
template<typename number>
struct __align__(16) mutual_trap {
	//number x, y, z;
	number stiff;
	number r0; 
	int p_ind;
};*/

template<typename number>
struct __align__(16) GPU_quat {
	number x, y, z, w;
};

/**
 * @brief Used when use_edge = true. It stores information associated to a single bond.
 */
typedef struct __align__(16) edge_bond {
	int from;
	int to;
	int n_from;
	int n_to;
} edge_bond;

/**
 * @brief Static class. It stores many utility functions used by CUDA classes. It could probably be turned into a namespace...
 */
class GpuUtils {
protected:
	static size_t _allocated_dev_mem;

public:
	template<typename number, typename number4>
	static number sum_4th_comp(number4 *v, int N) {
		number res = 0;
		for(int i = 0; i < N; i++) res += v[i].w;
		return res;
	}

	static float int_as_float(const int a) {
		union {
			int a;
			float b;
		} u;

		u.a = a;
		return u.b;
	}

	static int float_as_int(const float a) {
		union {
			float a;
			int b;
		} u;

		u.a = a;
		return u.b;
	}

#ifndef OLD_ARCH
	template<typename T> static void print_device_array(T *, int);
	template<typename T> static void check_device_thresold(T *, int, int);
#endif

	static size_t get_allocated_mem() { return _allocated_dev_mem; }
	static double get_allocated_mem_mb() { return get_allocated_mem() / 1048576.; }

	template<typename T>
	static cudaError_t LR_cudaMalloc(T **devPtr, size_t size);

};

template<typename T>
cudaError_t GpuUtils::LR_cudaMalloc(T **devPtr, size_t size) {
	GpuUtils::_allocated_dev_mem += size;
	return cudaMalloc(devPtr, size);
}

#endif /* GPUUTILS_H_ */
