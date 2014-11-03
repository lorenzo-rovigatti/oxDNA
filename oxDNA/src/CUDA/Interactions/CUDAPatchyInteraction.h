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

#define CUDA_MAX_PATCHES 4

/**
 * @brief CUDA implementation of the {@link PatchyInteraction patchy interaction}.
 */
template<typename number, typename number4>
class CUDAPatchyInteraction: public CUDABaseInteraction<number, number4>, public PatchyInteraction<number> {
public:
	CUDAPatchyInteraction();
	virtual ~CUDAPatchyInteraction();

	void get_settings(input_file &inp);
	void cuda_init(number box_side, int N);
	number get_cuda_rcut() { return this->get_rcut(); }

	void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds);
};

#endif /* CUDAPATCHYINTERACTION_H_ */
