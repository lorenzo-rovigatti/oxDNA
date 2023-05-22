/*
 * CUDAPolymerSwapInteraction.h
 *
 *  Created on: 24/mar/2020
 *      Author: lorenzo
 */

#ifndef CUDAPOLYMERSWAPINTERACTION_H_
#define CUDAPOLYMERSWAPINTERACTION_H_

#include "CUDA/Interactions/CUDABaseInteraction.h"
#include "PolymerSwapInteraction.h"

/**
 * @brief CUDA implementation of the {@link PolymerSwapInteraction interaction}.
 */
class CUDAPolymerSwapInteraction: public CUDABaseInteraction, public PolymerSwapInteraction {
private:
	c_number4 *_d_three_body_forces = nullptr;
	int *_d_bonded_neighs = nullptr;

public:
	CUDAPolymerSwapInteraction();
	virtual ~CUDAPolymerSwapInteraction();

	void get_settings(input_file &inp);
	void cuda_init(int N);
	c_number get_cuda_rcut() {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box);
};

extern "C" BaseInteraction *make_CUDAPolymerSwapInteraction() {
	return new CUDAPolymerSwapInteraction();
}

#endif /* CUDAPOLYMERSWAPINTERACTION_H_ */
