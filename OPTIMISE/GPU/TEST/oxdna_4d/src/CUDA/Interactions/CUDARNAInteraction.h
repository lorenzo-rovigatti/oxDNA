/*
 * CUDARNAInteraction.h
 *
 *  Created on: 22/jul/2014
 *      Author: petr
 */

#ifndef CUDARNAINTERACTION_H_
#define CUDARNAINTERACTION_H_

#include "CUDABaseInteraction.h"
#include "../../Interactions/RNAInteraction.h"

/**
 * @brief CUDA implementation of the oxRNA model, as provided by RNAInteraction.
 */

class CUDARNAInteraction: public CUDABaseInteraction, public RNAInteraction {
public:
	bool _grooving;
	bool _use_debye_huckel;
	bool _mismatch_repulsion;
	// copied from DNA2Interaction.h (CPU) (and change c_number -> float), the least bad way of doing things
	float _salt_concentration;
	bool _debye_huckel_half_charged_ends;
	float _debye_huckel_prefactor;
	float _debye_huckel_lambdafactor;
	float _RNA_HYDR_MIS;

	//the following values are calculated
	float _debye_huckel_RC; // this is the maximum interaction distance between backbones to interact with DH
	float _debye_huckel_RHIGH; // distance after which the potential is replaced by a quadratic cut-off
	float _debye_huckel_B; // prefactor of the quadratic cut-off
	float _minus_kappa;
	// End copy from DNA2Interaction.h
	int *_d_is_strand_end = nullptr;

	CUDARNAInteraction();
	virtual ~CUDARNAInteraction();

	void get_settings(input_file &inp) override;
	void cuda_init(int N) override;
	c_number get_cuda_rcut() {
		return this->get_rcut();
	}

	void _on_T_update() override;
	void _init_strand_ends(LR_bonds *d_bonds);

	void compute_forces(CUDABaseList*lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox*d_box);
	void _hb_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, float *hb_energies, int n_threads, bool *region_is_nearhb, CUDA_kernel_cfg hb_kernel_cfg, CUDABox*d_box);
	void _near_hb_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, bool *nearly_bonded_array, int n_threads, bool *region_is_nearhb, CUDA_kernel_cfg hb_kernel_cfg, CUDABox*d_box);
	void _dist_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, c_number *op_dists, int n_threads, CUDA_kernel_cfg dist_kernel_cfg, CUDABox*d_box);
};

#endif /* CUDARNAINTERACTION_H_ */
