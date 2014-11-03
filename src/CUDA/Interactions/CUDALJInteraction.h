/*
 * CUDALJInteraction.h
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#ifndef CUDALJINTERACTION_H_
#define CUDALJINTERACTION_H_

#include "CUDABaseInteraction.h"
#include "../../Interactions/LJInteraction.h"

/**
 * @brief CUDA implementation of the {@link LJInteraction Lennard-Jones interaction}.
 */
template<typename number, typename number4>
class CUDALJInteraction: public CUDABaseInteraction<number, number4>, public LJInteraction<number> {
public:
	CUDALJInteraction();
	virtual ~CUDALJInteraction();

	void get_settings(input_file &inp);
	void cuda_init(number box_side, int N);
	number get_cuda_rcut() { return this->get_rcut(); }

	void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds);
};

#endif /* CUDALJINTERACTION_H_ */
