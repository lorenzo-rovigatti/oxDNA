/*
 * cuda_hip_compat.h
 *
 * Maps the small set of CUDA runtime/texture symbols used by the oxDNA GPU
 * backend onto their HIP equivalents, so the .cu/.cuh sources compile with
 * hipcc without being modified. Active only when OXDNA_HIP is defined (set by
 * the HIP build path in CMake). The shim headers cuda.h / cuda_runtime.h /
 * cuda_runtime_api.h / vector_functions.h in this directory pull this in when
 * the existing sources #include the corresponding CUDA headers.
 */

#ifndef OXDNA_CUDA_HIP_COMPAT_H
#define OXDNA_CUDA_HIP_COMPAT_H

#ifdef OXDNA_HIP

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <hip/hip_runtime.h>

// status / types
#define cudaError_t                 hipError_t
#define cudaSuccess                 hipSuccess
#define cudaDeviceProp              hipDeviceProp_t

// memory management
#define cudaMalloc                  hipMalloc
#define cudaFree                    hipFree
#define cudaMemcpy                  hipMemcpy
#define cudaMemset                  hipMemset
#define cudaMemcpyHostToDevice      hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost      hipMemcpyDeviceToHost
#define cudaMemcpyDeviceToDevice    hipMemcpyDeviceToDevice

// HIP_SYMBOL(X) is X on AMD; wrapping keeps the NVIDIA path valid too.
#define cudaMemcpyToSymbol(symbol, ...) hipMemcpyToSymbol(HIP_SYMBOL(symbol), __VA_ARGS__)

// pinned host memory (oxDNA uses the nonstandard 3-arg cudaMallocHost form)
#define cudaMallocHost(ptr, ...)    hipHostMalloc((void **)(ptr), __VA_ARGS__)
#define cudaFreeHost                hipHostFree
#define cudaHostAllocDefault        hipHostMallocDefault

// device management
#define cudaSetDevice               hipSetDevice
#define cudaGetDevice               hipGetDevice
#define cudaGetDeviceCount          hipGetDeviceCount
#define cudaGetDeviceProperties     hipGetDeviceProperties
#define cudaDeviceSynchronize       hipDeviceSynchronize
#define cudaGetLastError            hipGetLastError
#define cudaGetErrorString          hipGetErrorString
#define cudaFuncCachePreferL1       hipFuncCachePreferL1

// AMD devices do not have a reconfigurable L1/shared-memory cache and
// hipDeviceSetCacheConfig returns hipErrorNotSupported on them. CUDA callers
// treat this as a required initialisation step, so preserve the CUDA call's
// intent as a successful no-op under HIP.
static inline hipError_t oxdna_hipDeviceSetCacheConfig(hipFuncCache_t) {
	return hipSuccess;
}
#define cudaDeviceSetCacheConfig    oxdna_hipDeviceSetCacheConfig

// texture objects (single use: the Verlet cell-counter lookup)
#define cudaTextureObject_t         hipTextureObject_t
#define cudaChannelFormatDesc       hipChannelFormatDesc
#define cudaResourceDesc            hipResourceDesc
#define cudaTextureDesc             hipTextureDesc
#define cudaCreateChannelDesc       hipCreateChannelDesc
#define cudaChannelFormatKindSigned hipChannelFormatKindSigned
#define cudaResourceTypeLinear      hipResourceTypeLinear
#define cudaReadModeElementType     hipReadModeElementType
#define cudaCreateTextureObject     hipCreateTextureObject
#define cudaDestroyTextureObject    hipDestroyTextureObject

// The only symbol the trimmed NVIDIA helper_cuda.h actually exports; under HIP
// that header is inert (its blocks are guarded by CUDA-only macros), so provide
// it here.
static inline void __oxdna_getLastCudaError(const char *msg, const char *file, int line) {
	hipError_t err = hipGetLastError();
	if(err != hipSuccess) {
		printf("HIP error at %s:%d : %s : %s\n", file, line, msg, hipGetErrorString(err));
		exit(EXIT_FAILURE);
	}
}
#ifndef getLastCudaError
#define getLastCudaError(msg) __oxdna_getLastCudaError((msg), __FILE__, __LINE__)
#endif

// oxDNA calls sincos() with float arguments; HIP provides sincos(double,...) and
// sincosf(float,...) but not sincos(float,...). Supply the float overload.
__device__ static inline void sincos(float x, float *s, float *c) { sincosf(x, s, c); }

// oxDNA writes pointer parameters as `type __restrict__ *ptr` (the qualifier
// before the star), which nvcc tolerates but clang rejects. __restrict__ is only
// an aliasing hint, so drop it for HIP. (Performance follow-up: restore a
// well-formed `* __restrict__` form to keep the hint.)
#define __restrict__

#endif // OXDNA_HIP
#endif // OXDNA_CUDA_HIP_COMPAT_H
