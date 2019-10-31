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

class CUDALJInteraction: public CUDABaseInteraction, public LJInteraction {
public:
	CUDALJInteraction();
	virtual ~CUDALJInteraction();

	void get_settings(input_file &inp);
	void cuda_init(c_number box_side, int N);
	c_number get_cuda_rcut() {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList*lists, tmpnmbr *d_poss, GPU_quat *d_orientations, tmpnmbr *d_forces, tmpnmbr *d_torques, LR_bonds *d_bonds, CUDABox*d_box);
};

#endif /* CUDALJINTERACTION_H_ */
