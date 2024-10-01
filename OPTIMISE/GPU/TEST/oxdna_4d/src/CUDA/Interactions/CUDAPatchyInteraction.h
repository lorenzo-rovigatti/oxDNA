/*
 * CUDAPatchyInteraction.h
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#ifndef CUDAPATCHYINTERACTION_H_
#define CUDAPATCHYINTERACTION_H_

#include "CUDABaseInteraction.h"
#include "../../Interactions/PatchyInteraction.h"

#define CUDA_MAX_PATCHES 5

/**
 * @brief CUDA implementation of the {@link PatchyInteraction patchy interaction}.
 */

class CUDAPatchyInteraction: public CUDABaseInteraction, public PatchyInteraction {
public:
	CUDAPatchyInteraction();
	virtual ~CUDAPatchyInteraction();

	void get_settings(input_file &inp) override;
	void cuda_init(int N) override;
	c_number get_cuda_rcut() {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList*lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox*d_box);
};

#endif /* CUDAPATCHYINTERACTION_H_ */
