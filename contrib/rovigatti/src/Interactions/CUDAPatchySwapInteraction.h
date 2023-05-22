/*
 * CUDAPatchySwapInteraction.h
 *
 *  Created on: 15/jul/2020
 *      Author: lorenzo
 */

#ifndef CUDAPATCHYSWAPINTERACTION_H_
#define CUDAPATCHYSWAPINTERACTION_H_

#include "CUDA/Interactions/CUDABaseInteraction.h"

#include "PatchySwapInteraction.h"

struct CUDA_FS_bonding_pattern;
struct swap_event;

/**
 * @brief CUDA implementation of the {@link PatchySwapInteraction}.
 */

class CUDAPatchySwapInteraction: public CUDABaseInteraction, public PatchySwapInteraction {
protected:
	c_number4 *_d_three_body_forces, *_d_three_body_torques;
	llint _step;
public:
	static const int MAX_PATCHES = 5;
	static const int MAX_SPECIES = 5;
	static const int MAX_NEIGHS = 20;

	CUDAPatchySwapInteraction();
	virtual ~CUDAPatchySwapInteraction();

	void get_settings(input_file &inp);
	void cuda_init(int N);
	c_number get_cuda_rcut() {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box);
};

extern "C" BaseInteraction *make_CUDAPatchySwapInteraction() {
	return new CUDAPatchySwapInteraction();
}

#endif /* CUDAPATCHYSWAPINTERACTION_H_ */
