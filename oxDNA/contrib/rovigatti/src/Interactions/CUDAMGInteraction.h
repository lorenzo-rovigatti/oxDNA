/*
 * CUDACPMixtureInteraction.h
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#ifndef CUDAMGINTERACTION_H_
#define CUDAMGINTERACTION_H_

#include "CUDA/Interactions/CUDABaseInteraction.h"
#include "MGInteraction.h"

/**
 * @brief CUDA implementation of the {@link MGInteraction MicroGel interaction}.
 */
template<typename number, typename number4>
class CUDAMGInteraction: public CUDABaseInteraction<number, number4>, public MGInteraction<number> {
private:
	int *_d_bonded_neighs;
public:
	CUDAMGInteraction();
	virtual ~CUDAMGInteraction();

	void get_settings(input_file &inp);
	void cuda_init(number box_side, int N);
	number get_cuda_rcut() { return this->get_rcut(); }

	void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box);
};

extern "C" IBaseInteraction<float> *make_CUDAMGInteraction_float() { return new CUDAMGInteraction<float, float4>(); }
extern "C" IBaseInteraction<double> *make_CUDAMGInteraction_double() { return new CUDAMGInteraction<double, LR_double4>(); }

#endif /* CUDAMGINTERACTION_H_ */
