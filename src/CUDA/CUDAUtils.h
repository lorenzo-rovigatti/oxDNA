/**
 * @file    CUDAUtils.h
 * @date    24/set/2010
 * @author  lorenzo
 *
 *
 */

#ifndef GPUUTILS_H_
#define GPUUTILS_H_

#include "../defs.h"
#include "cuda_defs.h"

#include "cuda_utils/helper_cuda.h"

/**
* @brief Static class. It stores many utility functions used by CUDA classes. It could probably be turned into a namespace...
*/
class GpuUtils {
protected:
	static size_t _allocated_dev_mem;

public:
	template<typename T> static void print_device_array(T *, int);
	template<typename T> static void check_device_thresold(T *, int, int);

	static c_number sum_4th_comp(c_number4 *v, int N) {
		c_number res = 0;
		for(int i = 0; i < N; i++)
			res += v[i].w;
		return res;
	}

	static c_number4 sum_c_number4_on_GPU(c_number4 *dv, int N);

	static double sum_c_number4_to_double_on_GPU(c_number4 *dv, int N);

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

	static size_t get_allocated_mem() {
		return _allocated_dev_mem;
	}
	static double get_allocated_mem_mb() {
		return get_allocated_mem() / 1048576.;
	}

	template<typename T>
	static cudaError_t LR_cudaMalloc(T **devPtr, size_t size);

	static void reset_allocated_mem() {
		_allocated_dev_mem = 0;
	}

	template<typename T>
	static void init_texture_object(cudaTextureObject_t *obj, cudaChannelFormatDesc format, T *dev_ptr, size_t size) {
		cudaResourceDesc res_desc_eps;
		memset(&res_desc_eps, 0, sizeof(res_desc_eps));
		res_desc_eps.resType = cudaResourceTypeLinear;
		res_desc_eps.res.linear.devPtr = dev_ptr;
		res_desc_eps.res.linear.desc = format;
		res_desc_eps.res.linear.sizeInBytes = size * sizeof(T);

		cudaTextureDesc tex_desc_eps;
		memset(&tex_desc_eps, 0, sizeof(tex_desc_eps));
		tex_desc_eps.readMode = cudaReadModeElementType;

		CUDA_SAFE_CALL(cudaCreateTextureObject(obj, &res_desc_eps, &tex_desc_eps, NULL));
	}
};

template<typename T>
cudaError_t GpuUtils::LR_cudaMalloc(T **devPtr, size_t size) {
	GpuUtils::_allocated_dev_mem += size;
	return cudaMalloc((void **) devPtr, size);
}

#endif /* GPUUTILS_H_ */
