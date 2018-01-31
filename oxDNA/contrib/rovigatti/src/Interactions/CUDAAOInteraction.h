/*
 * CUDAAOInteraction.h
 *
 *  Created on: 25/oct/2017
 *      Author: lorenzo
 */

#ifndef CUDAAOINTERACTION_H_
#define CUDAAOINTERACTION_H_

#include "CUDA/Interactions/CUDABaseInteraction.h"
#include "AOInteraction.h"

/**
 * @brief CUDA implementation of the {@link AOInteraction Asakura-Oosawa interaction}.
 */
template<typename number, typename number4>
class CUDAAOInteraction: public CUDABaseInteraction<number, number4>, public AOInteraction<number> {
public:
	CUDAAOInteraction();
	virtual ~CUDAAOInteraction();

	void get_settings(input_file &inp);
	void cuda_init(number box_side, int N);
	number get_cuda_rcut() { return this->get_rcut(); }

	void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box);
};

extern "C" IBaseInteraction<float> *make_CUDAAOInteraction_float() { return new CUDAAOInteraction<float, float4>(); }
extern "C" IBaseInteraction<double> *make_CUDAAOInteraction_double() { return new CUDAAOInteraction<double, LR_double4>(); }

#endif /* CUDAAOINTERACTION_H_ */
