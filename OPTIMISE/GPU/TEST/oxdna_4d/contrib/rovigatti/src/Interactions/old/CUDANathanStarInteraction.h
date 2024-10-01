/*
 * CUDANathanStarInteraction.h
 *
 *  Created on: 10/apr/2015
 *      Author: lorenzo
 */

#ifndef CUDANATHANSTARINTERACTION_H_
#define CUDANATHANSTARINTERACTION_H_

#include "CUDA/Interactions/CUDABaseInteraction.h"

#include "NathanStarInteraction.h"

/**
 * @brief CUDA implementation of the {@link NathanStarInteraction}.
 */

class CUDANathanStarInteraction: public CUDABaseInteraction, public NathanStarInteraction {
protected:
	float *_d_patchy_star;

	void _setup_cuda_interp();
public:
	CUDANathanStarInteraction();
	virtual ~CUDANathanStarInteraction();

	void get_settings(input_file &inp);
	void cuda_init(int N);
	c_number get_cuda_rcut() {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box);
};

extern "C" BaseInteraction *make_CUDANathanStarInteraction() {
	return new CUDANathanStarInteraction();
}

#endif /* CUDANATHANSTARINTERACTION_H_ */
