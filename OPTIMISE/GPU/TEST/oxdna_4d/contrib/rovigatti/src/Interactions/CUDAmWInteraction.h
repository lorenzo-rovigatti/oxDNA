/*
 * CUDAmWInteraction.h
 *
 *  Created on: 10/apr/2015
 *      Author: lorenzo
 */

#ifndef CUDAMWINTERACTION_H_
#define CUDAMWINTERACTION_H_

#include "CUDA/Interactions/CUDABaseInteraction.h"

#include "mWInteraction.h"

/**
 * @brief CUDA implementation of the {@link mWInteraction}.
 */

class CUDAmWInteraction: public CUDABaseInteraction, public mWInteraction {
protected:
	c_number4 *_d_forces_3body;
public:
	CUDAmWInteraction();
	virtual ~CUDAmWInteraction();

	void get_settings(input_file &inp);
	void cuda_init(int N);
	c_number get_cuda_rcut() {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box);
};

extern "C" BaseInteraction *make_CUDAmWInteraction() {
	return new CUDAmWInteraction();
}

#endif /* CUDAMWINTERACTION_H_ */
