/*
 * CUDADetailedPatchySwapInteraction.h
 *
 *  Created on: 14/may/2021
 *      Author: lorenzo
 */

#ifndef CUDADETAILEDPATCHYSWAPINTERACTION_H_
#define CUDADETAILEDPATCHYSWAPINTERACTION_H_

#include "CUDA/Interactions/CUDABaseInteraction.h"

#include "DetailedPatchySwapInteraction.h"

struct CUDA_FS_bonding_pattern;
struct swap_event;

/**
 * @brief CUDA implementation of the {@link PatchySwapInteraction}.
 */

class CUDADetailedPatchySwapInteraction: public CUDABaseInteraction, public DetailedPatchySwapInteraction {
protected:
	c_number4 *_d_three_body_forces = nullptr;
	c_number4 *_d_three_body_torques = nullptr;

	float *_d_patchy_eps = nullptr;
	float4 *_d_base_patches = nullptr;

	cudaTextureObject_t _tex_patchy_eps = 0;
	cudaTextureObject_t _tex_base_patches = 0;

	llint _step;
public:
	static const int MAX_PATCHES = 5;
	static const int MAX_SPECIES = 60;
	static const int MAX_NEIGHS = 20;

	CUDADetailedPatchySwapInteraction();
	virtual ~CUDADetailedPatchySwapInteraction();

	void get_settings(input_file &inp);
	void cuda_init(int N);
	c_number get_cuda_rcut() {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box);
};

extern "C" BaseInteraction *make_CUDADetailedPatchySwapInteraction() {
	return new CUDADetailedPatchySwapInteraction();
}

#endif /* CUDAPATCHYSWAPINTERACTION_H_ */
