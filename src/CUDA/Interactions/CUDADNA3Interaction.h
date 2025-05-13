/*
 * CUDADNA3Interaction.h
 *
 *  Created on: 13/may/25
 *      Author: lorenzo
 */

#ifndef CUDADNA3INTERACTION_H_
#define CUDADNA3INTERACTION_H_

#include "CUDABaseInteraction.h"
#include "../../Interactions/DNAInteraction.h"

/**
 * @brief CUDA implementation of the oxDNA3 model, as provided by DNA3Interaction.
 */
class CUDADNA3Interaction: public CUDABaseInteraction, public DNAInteraction {
public:
	enum {
		DEBYE_HUCKEL = 7
	};
	CUDADNA3Interaction();
	virtual ~CUDADNA3Interaction();

	bool _use_debye_huckel;
	bool _use_oxDNA2_coaxial_stacking;
	bool _use_oxDNA2_FENE;
	// copied from DNA2Interaction.h (CPU) (and change c_number -> float), the least bad way of doing things
	float _salt_concentration;
	bool _debye_huckel_half_charged_ends;
	float _debye_huckel_prefactor;
	float _debye_huckel_lambdafactor;

	//the following values are calculated
	float _debye_huckel_RC; // this is the maximum interaction distance between backbones to interact with DH
	float _debye_huckel_RHIGH; // distance after which the potential is replaced by a quadratic cut-off
	float _debye_huckel_B; // prefactor of the quadratic cut-off
	float _minus_kappa;
	// End copy from DNA2Interaction.h
	int *_d_is_strand_end = nullptr;

	void get_settings(input_file &inp) override;
	void cuda_init(int N) override;
	c_number get_cuda_rcut() {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList*lists, c_number4 *d_poss, GPU_quat *d_qorientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox*d_box);

protected:
	void _on_T_update() override;
	void _init_strand_ends(LR_bonds *d_bonds);
};

#endif /* CUDADNA3INTERACTION_H_ */
