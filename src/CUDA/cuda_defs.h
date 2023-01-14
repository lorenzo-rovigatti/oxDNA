/*
 * cuda_defs.h
 *
 *  Created on: 25 ott 2019
 *      Author: lorenzo
 */

#ifndef SRC_CUDA_CUDA_DEFS_H_
#define SRC_CUDA_CUDA_DEFS_H_

/// CUDA_SAFE_CALL replacement for backwards compatibility (CUDA < 5.0)
#define CUDA_SAFE_CALL(call)                                  \
  do {                                                        \
    cudaError_t err = call;                                   \
    if (err != cudaSuccess) {                                 \
      printf("CUDA error at %s %d: %s\n", __FILE__, __LINE__, \
             cudaGetErrorString(err));                        \
      exit(EXIT_FAILURE);                                     \
    }                                                         \
  } while (0)
/// CUT_CHECK_ERROR replacement for backwards compatibility (CUDA < 5.0)
#define CUT_CHECK_ERROR(x) getLastCudaError(x);

/// threads per block
#define TINBLOCK (blockDim.x*blockDim.y)
/// c_number of blocks
#define NBLOCKS (gridDim.x*gridDim.y)
/// c_number of threads
#define NTHREADS (NBLOCKS * TINBLOCK)

/// thread id relative to its block
#define TID (blockDim.x*threadIdx.y + threadIdx.x)
/// block id
#define BID (gridDim.x*blockIdx.y + blockIdx.x)
/// thread id
#define IND (TINBLOCK * BID + TID)

#define CUDA_LRACOS(x) (((x) >= (c_number)1) ? (c_number) 0 : ((x) <= (c_number)-1) ? (c_number) PI : acosf(x))
#define CUDA_DOT(a, b) (a.x*b.x + a.y*b.y + a.z*b.z)

#define COPY_ARRAY_TO_CONSTANT(dest, src, size) {\
		float *val = new float[(size)];\
		for(int i = 0; i < (size); i++) val[i] = (float) ((src)[i]);\
		CUDA_SAFE_CALL(cudaMemcpyToSymbol((dest), val, (size)*sizeof(float)));\
		delete[] val; }

#define COPY_NUMBER_TO_FLOAT(dest, src) {\
		float tmp = src;\
		CUDA_SAFE_CALL(cudaMemcpyToSymbol((dest), &tmp, sizeof(float)));\
		}

#include <cuda.h>
#include <cuda_runtime_api.h>

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
typedef struct
__align__(16) {
	double x, y, z;
	float w;
} LR_double4;

#ifdef CUDA_DOUBLE_PRECISION
using c_number4 = LR_double4;
using c_number = double;
using GPU_quat = double4;
#else
using c_number4 = float4;
using c_number = float;
using GPU_quat = float4;
#endif

/**
 * @brief It keeps track of neighbours along 3" and 5" directions.
 */
typedef struct
__align__(8) {
	int n3, n5;
} LR_bonds;

/**
 * @brief Used when use_edge = true. It stores information associated to a single bond.
 */
typedef struct __align__(8) edge_bond {
	int from;
	int to;
} edge_bond;

/**
 * @brief Used to store the stress tensor on GPUs
 */
typedef struct __align__(16) CUDAStressTensor {
	c_number e[6];

	__device__ __host__ CUDAStressTensor() : e{0} {

	}

	__device__ __host__ CUDAStressTensor(c_number e0, c_number e1, c_number e2, c_number e3, c_number e4, c_number e5) :
		e{e0, e1, e2, e3, e4, e5} {

	}

	__device__ __host__ inline CUDAStressTensor operator+(const CUDAStressTensor &other) const {
		return CUDAStressTensor(
				e[0] + other.e[0],
				e[1] + other.e[1],
				e[2] + other.e[2],
				e[3] + other.e[3],
				e[4] + other.e[4],
				e[5] + other.e[5]);
	}

	__device__ __host__ inline void operator+=(const CUDAStressTensor &other) {
		e[0] += other.e[0];
		e[1] += other.e[1];
		e[2] += other.e[2];
		e[3] += other.e[3];
		e[4] += other.e[4];
		e[5] += other.e[5];
	}

	__device__ __host__ inline void operator*=(const c_number other) {
		e[0] *= other;
		e[1] *= other;
		e[2] *= other;
		e[3] *= other;
		e[4] *= other;
		e[5] *= other;
	}

	__device__ __host__ inline void operator/=(const c_number other) {
		e[0] /= other;
		e[1] /= other;
		e[2] /= other;
		e[3] /= other;
		e[4] /= other;
		e[5] /= other;
	}

	__host__ StressTensor as_StressTensor() {
		return StressTensor({e[0], e[1], e[2], e[3], e[4], e[5]});
	}
} CUDAStressTensor;

#endif /* SRC_CUDA_CUDA_DEFS_H_ */
