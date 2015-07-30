/*
 * CUDANathanStarInteraction.h
 *
 *  Created on: 10/apr/2015
 *      Author: lorenzo
 */

#ifndef CUDANATHANSTARINTERACTION_H_
#define CUDANATHANSTARINTERACTION_H_

#include "CUDABaseInteraction.h"

#include "NathanStarInteraction.h"

/**
 * @brief CUDA implementation of the {@link NathanStarInteraction}.
 */
template<typename number, typename number4>
class CUDANathanStarInteraction: public CUDABaseInteraction<number, number4>, public NathanStarInteraction<number> {
protected:
	float *_d_patchy_star;

	void _setup_cuda_interp();
public:
	CUDANathanStarInteraction();
	virtual ~CUDANathanStarInteraction();

	void get_settings(input_file &inp);
	void cuda_init(number box_side, int N);
	number get_cuda_rcut() { return this->get_rcut(); }

	void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds);
};

extern "C" IBaseInteraction<float> *make_CUDANathanStarInteraction_float() { return new CUDANathanStarInteraction<float, float4>(); }
extern "C" IBaseInteraction<double> *make_CUDANathanStarInteraction_double() { return new CUDANathanStarInteraction<double, LR_double4>(); }

#endif /* CUDANATHANSTARINTERACTION_H_ */
