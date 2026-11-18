/*
 * CUDACGNucleicAcidsInteraction.h
 *
 *  Created on: 24/mar/2020
 *      Author: lorenzo
 */

#ifndef CUDACGNUCLEICACIDSINTERACTION_H_
#define CUDACGNUCLEICACIDSINTERACTION_H_

#include "CUDA/Interactions/CUDABaseInteraction.h"
#include "ANNaMoInteraction.h"

/**
 * @brief CUDA implementation of the {@link CGNucleicAcidsInteraction interaction}.
 */
class CUDACGNucleicAcidsInteraction: public CUDABaseInteraction, public CGNucleicAcidsInteraction {
private:
	/// copies the parameters of the capped FENE potential (max_backbone_force) to the GPU
	void _copy_mbf_to_device();

	c_number4 *_d_three_body_forces = nullptr;
	c_number4 *_d_three_body_torques = nullptr;
	float *_d_3b_epsilon = nullptr;
	int *_d_crossovers = nullptr;

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

/// the entry point used by the CUDA factory for interaction_type = ANNaMoInteraction
extern "C" BaseInteraction *make_CUDAANNaMoInteraction() {
	return new CUDACGNucleicAcidsInteraction();
}

#endif /* CUDACGNUCLEICACIDSINTERACTION_H_ */
