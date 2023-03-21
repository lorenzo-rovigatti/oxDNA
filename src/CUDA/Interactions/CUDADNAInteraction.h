/*
 * CUDADNAInteraction.h
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#ifndef CUDADNAINTERACTION_H_
#define CUDADNAINTERACTION_H_

#include "CUDABaseInteraction.h"
#include "../../Interactions/DNAInteraction.h"

/**
 * @brief CUDA implementation of the oxDNA model, as provided by DNAInteraction.
 */

class CUDADNAInteraction: public CUDABaseInteraction, public DNAInteraction {
public:
	enum {
		DEBYE_HUCKEL = 7
	};
	CUDADNAInteraction();
	virtual ~CUDADNAInteraction();

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
	void _hb_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, float *hb_energies, int n_threads, bool *region_is_nearhb, CUDA_kernel_cfg hb_kernel_cfg, CUDABox*d_box);
	void _near_hb_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, bool *nearly_bonded_array, int n_threads, bool *region_is_nearhb, CUDA_kernel_cfg hb_kernel_cfg, CUDABox*d_box);
	void _dist_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, c_number *op_dists, int n_threads, CUDA_kernel_cfg dist_kernel_cfg, CUDABox*d_box);

protected:
	void _on_T_update() override;
	void _init_strand_ends(LR_bonds *d_bonds);
};

#endif /* CUDADNAINTERACTION_H_ */
