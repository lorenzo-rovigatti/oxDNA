/*
 * CUDADetailedPolymerSwapInteraction.h
 *
 *  Created on: 17/mar/2022
 *      Author: lorenzo
 */

#ifndef CUDADETAILEDPOLYMERSWAPINTERACTION_H_
#define CUDADETAILEDPOLYMERSWAPINTERACTION_H_

#include "CUDA/Interactions/CUDABaseInteraction.h"
#include "DetailedPolymerSwapInteraction.h"

/**
 * @brief CUDA implementation of the {@link DetailedPolymerSwapInteraction interaction}.
 */

class CUDADetailedPolymerSwapInteraction: public CUDABaseInteraction, public DetailedPolymerSwapInteraction {
private:
	c_number4 *_d_three_body_forces;
	int *_d_bonded_neighs;
	float *_d_3b_epsilon = nullptr;

	cudaTextureObject_t _tex_eps = 0;

public:
	CUDADetailedPolymerSwapInteraction();
	virtual ~CUDADetailedPolymerSwapInteraction();

	void get_settings(input_file &inp);
	void cuda_init(int N);
	c_number get_cuda_rcut() {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box);
};

extern "C" BaseInteraction *make_CUDADetailedPolymerSwapInteraction() {
	return new CUDADetailedPolymerSwapInteraction();
}

#endif /* CUDADETAILEDPOLYMERSWAPINTERACTION_H_ */
