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
template<typename number, typename number4>
class CUDADNAInteraction: public CUDABaseInteraction<number, number4>, public DNAInteraction<number> {
public:
        enum {
		DEBYE_HUCKEL = 7
	};
	CUDADNAInteraction();
	virtual ~CUDADNAInteraction();

	bool _use_debye_huckel;
	bool _use_oxDNA2_coaxial_stacking;
	bool _use_oxDNA2_FENE;
	// copied from DNA2Interaction.h (CPU) (and change number -> float), the least bad way of doing things
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

	void get_settings(input_file &inp);
	void cuda_init(number box_side, int N);
	number get_cuda_rcut() { return this->get_rcut(); }

	void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_qorientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds);
	void _hb_op_precalc(number4 *poss, GPU_quat<number> *orientations, int *op_pairs1, int *op_pairs2, float *hb_energies, int n_threads, bool *region_is_nearhb, CUDA_kernel_cfg hb_kernel_cfg);
	void _near_hb_op_precalc(number4 *poss, GPU_quat<number> *orientations, int *op_pairs1, int *op_pairs2, bool *nearly_bonded_array, int n_threads, bool *region_is_nearhb, CUDA_kernel_cfg hb_kernel_cfg);
	void _dist_op_precalc(number4 *poss, GPU_quat<number> *orientations, int *op_pairs1, int *op_pairs2, number *op_dists, int n_threads, CUDA_kernel_cfg dist_kernel_cfg);
};

#endif /* CUDADNAINTERACTION_H_ */
