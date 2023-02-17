/*
 * CUDAFSInteraction.h
 *
 *  Created on: 10/apr/2015
 *      Author: lorenzo
 */

#ifndef CUDAFSINTERACTION_H_
#define CUDAFSINTERACTION_H_

#include "CUDA/Interactions/CUDABaseInteraction.h"

#include "FSInteraction.h"

struct CUDA_FS_bonding_pattern;
struct swap_event;

/**
 * @brief CUDA implementation of the {@link FSInteraction}.
 */

class CUDAFSInteraction: public CUDABaseInteraction, public FSInteraction {
protected:
	c_number4 *_d_three_body_forces, *_d_three_body_torques;
	llint _step;

	// polymer stuff
	int *_d_bonded_neighs;
public:
	CUDAFSInteraction();
	virtual ~CUDAFSInteraction();

	void get_settings(input_file &inp);
	void cuda_init(int N);
	c_number get_cuda_rcut() {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box);
};

extern "C" BaseInteraction *make_CUDAFSInteraction() {
	return new CUDAFSInteraction();
}

#endif /* CUDAFSINTERACTION_H_ */
