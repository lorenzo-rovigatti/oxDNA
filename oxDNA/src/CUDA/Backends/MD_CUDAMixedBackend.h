/*
 * MD_CUDAMixedBackend.h
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#ifndef MD_CUDAMIXEDBACKEND_H_
#define MD_CUDAMIXEDBACKEND_H_

#include "MD_CUDABackend.h"

/**
 * @brief CUDA backend with mixed precision for MD simulations.
 *
 * This class is a regular MD backend written to be almost as fast as a MD_CUDABackend<float, float4>
 * and, at the same time, almost as reliable as MD_CUDABackend<double, LR_double4> when it comes
 * to numerical precision. This is probably the best class for production simulations.
 */
class CUDAMixedBackend: public MD_CUDABackend<float, float4> {
protected:
	LR_double4 *_d_possd;
	LR_double4 *_d_velsd, *_d_Lsd;
	GPU_quat<double> *_d_orientationsd;
	size_t _vec_sized, _orient_sized;

	void _float4_to_LR_double4(float4 *src, LR_double4 *dest);
	void _LR_double4_to_float4(LR_double4 *src, float4 *dest);
	void _quat_double_to_quat_float(GPU_quat<double> *src, GPU_quat<float> *dest);
	void _quat_float_to_quat_double(GPU_quat<float> *src, GPU_quat<double> *dest);

	void _init_CUDA_MD_symbols();

	virtual void _sort_particles();

	virtual void _first_step();
	virtual void _forces_second_step();

	virtual void _thermalize(llint curr_step);
	virtual void _gpu_to_host_particles();
	virtual void _host_particles_to_gpu();

public:
	CUDAMixedBackend();
	virtual ~CUDAMixedBackend();

	void init();
};

#endif /* MD_CUDAMIXEDBACKEND_H_ */
