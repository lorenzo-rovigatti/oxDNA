/*
 * CUDACGNucleicAcidsInteraction.h
 *
 *  Created on: 24/mar/2020
 *      Author: lorenzo
 */

#ifndef CUDACGNUCLEICACIDSINTERACTION_H_
#define CUDACGNUCLEICACIDSINTERACTION_H_

#include "CUDA/Interactions/CUDABaseInteraction.h"
#include "CGNucleicAcidsInteraction.h"

/**
 * @brief CUDA implementation of the {@link CGNucleicAcidsInteraction interaction}.
 */
class CUDACGNucleicAcidsInteraction: public CUDABaseInteraction, public CGNucleicAcidsInteraction {
private:
	c_number4 *_d_three_body_forces = nullptr;
	int *_d_bonded_neighs = nullptr;
	float *_d_3b_epsilon = nullptr;

	cudaTextureObject_t _tex_eps = 0;

public:
	CUDACGNucleicAcidsInteraction();
	virtual ~CUDACGNucleicAcidsInteraction();

	void get_settings(input_file &inp);
	void cuda_init(int N);
	c_number get_cuda_rcut() {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box);
};

extern "C" BaseInteraction *make_CUDACGNucleicAcidsInteraction() {
	return new CUDACGNucleicAcidsInteraction();
}

#endif /* CUDACGNUCLEICACIDSINTERACTION_H_ */
