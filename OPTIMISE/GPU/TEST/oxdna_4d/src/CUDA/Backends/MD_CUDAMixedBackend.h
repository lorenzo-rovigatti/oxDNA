/*
 * MD_CUDAMixedBackend.h
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#ifndef MD_CUDAMIXEDBACKEND_H_
#define MD_CUDAMIXEDBACKEND_H_

#include "MD_CUDABackend.h"

// this compilation unit is compiled only when the CUDA backend is compiled with CUDA_DOUBLE set to OFF
using GPU_quat_double = double4;

/**
 * @brief CUDA backend with mixed precision for MD simulations.
 *
 * This class is a regular MD backend written to be almost as fast as a MD_CUDABackend<float, float4>
 * and, at the same time, almost as reliable as MD_CUDABackend<double, LR_double4> when it comes
 * to numerical precision. This is probably the best class for production simulations.
 */
class CUDAMixedBackend: public MD_CUDABackend {
protected:
	LR_double4 *_d_possd;
	LR_double4 *_d_velsd, *_d_Lsd;
	GPU_quat_double *_d_orientationsd;
	size_t _vec_sized = 0;
	size_t _orient_sized = 0;

	void _float4_to_LR_double4(float4 *src, LR_double4 *dest);
	void _LR_double4_to_float4(LR_double4 *src, float4 *dest);
	void _quat_double_to_quat_float(GPU_quat_double *src, GPU_quat *dest);
	void _quat_float_to_quat_double(GPU_quat *src, GPU_quat_double *dest);

	void _init_CUDA_MD_symbols() override;

	virtual void _sort_particles();
	virtual void _rescale_positions(float4 new_Ls, float4 old_Ls);

	virtual void _first_step();
	virtual void _forces_second_step();

	void _thermalize() override;
	void _update_stress_tensor() override;

public:
	CUDAMixedBackend();
	virtual ~CUDAMixedBackend();

	void init();

	virtual void apply_simulation_data_changes();
	virtual void apply_changes_to_simulation_data();
};

#endif /* MD_CUDAMIXEDBACKEND_H_ */
