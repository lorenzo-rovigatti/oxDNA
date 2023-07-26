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

class CUDAAOInteraction: public CUDABaseInteraction, public AOInteraction {
public:
	CUDAAOInteraction();
	virtual ~CUDAAOInteraction();

	void get_settings(input_file &inp);
	void cuda_init(int N);
	c_number get_cuda_rcut() {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box);
};

extern "C" BaseInteraction *make_CUDAAOInteraction() {
	return new CUDAAOInteraction();
}

#endif /* CUDAAOINTERACTION_H_ */
