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

class CUDACPMixtureInteraction: public CUDABaseInteraction, public CPMixtureInteraction {
public:
	CUDACPMixtureInteraction();
	virtual ~CUDACPMixtureInteraction();

	void get_settings(input_file &inp);
	void cuda_init(int N);
	c_number get_cuda_rcut() {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box);
};

extern "C" BaseInteraction *make_CUDACPMixtureInteraction() {
	return new CUDACPMixtureInteraction();
}

#endif /* CUDACPMIXTUREINTERACTION_H_ */
