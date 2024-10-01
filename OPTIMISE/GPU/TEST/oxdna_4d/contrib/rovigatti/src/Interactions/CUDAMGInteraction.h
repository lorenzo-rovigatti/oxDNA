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

class CUDAMGInteraction: public CUDABaseInteraction, public MGInteraction {
private:
	int *_d_bonded_neighs;
public:
	CUDAMGInteraction();
	virtual ~CUDAMGInteraction();

	void get_settings(input_file &inp);
	void cuda_init(int N);
	c_number get_cuda_rcut() {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box);
};

extern "C" BaseInteraction *make_CUDAMGInteraction() {
	return new CUDAMGInteraction();
}

#endif /* CUDAMGINTERACTION_H_ */
