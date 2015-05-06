/*
 * CUDAmWInteraction.h
 *
 *  Created on: 10/apr/2015
 *      Author: lorenzo
 */

#ifndef CUDAMWINTERACTION_H_
#define CUDAMWINTERACTION_H_

#include "CUDABaseInteraction.h"

#include "mWInteraction.h"

/**
 * @brief CUDA implementation of the {@link mWInteraction}.
 */
template<typename number, typename number4>
class CUDAmWInteraction: public CUDABaseInteraction<number, number4>, public mWInteraction<number> {
protected:
	number4 *_d_forces_3body;
public:
	CUDAmWInteraction();
	virtual ~CUDAmWInteraction();

	void get_settings(input_file &inp);
	void cuda_init(number box_side, int N);
	number get_cuda_rcut() { return this->get_rcut(); }

	void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds);
};

extern "C" IBaseInteraction<float> *make_CUDAmWInteraction_float() { return new CUDAmWInteraction<float, float4>(); }
extern "C" IBaseInteraction<double> *make_CUDAmWInteraction_double() { return new CUDAmWInteraction<double, LR_double4>(); }

#endif /* CUDAMWINTERACTION_H_ */
