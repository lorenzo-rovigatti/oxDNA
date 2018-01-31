/*
 * CUDACPMixtureInteraction.h
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#ifndef CUDACPMIXTUREINTERACTION_H_
#define CUDACPMIXTUREINTERACTION_H_

#include "CUDA/Interactions/CUDABaseInteraction.h"
#include "CPMixtureInteraction.h"

/**
 * @brief CUDA implementation of the {@link CPMixtureInteraction Colloid-polymer mixture interaction}.
 */
template<typename number, typename number4>
class CUDACPMixtureInteraction: public CUDABaseInteraction<number, number4>, public CPMixtureInteraction<number> {
public:
	CUDACPMixtureInteraction();
	virtual ~CUDACPMixtureInteraction();

	void get_settings(input_file &inp);
	void cuda_init(number box_side, int N);
	number get_cuda_rcut() { return this->get_rcut(); }

	void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box);
};

extern "C" IBaseInteraction<float> *make_CUDACPMixtureInteraction_float() { return new CUDACPMixtureInteraction<float, float4>(); }
extern "C" IBaseInteraction<double> *make_CUDACPMixtureInteraction_double() { return new CUDACPMixtureInteraction<double, LR_double4>(); }

#endif /* CUDACPMIXTUREINTERACTION_H_ */
