/*
 * CUDADNA2ModInteraction.h
 *
 *  Created on: 5/feb/2018
 *      Author: Ferdinando 
 *      re-Author: Ferdinando, after CUDADNAInteraction.h by Lorenzo
 */

#ifndef CUDADNA2MODINTERACTION_H_
#define CUDADNA2MODINTERACTION_H_

#include "CUDA/Interactions/CUDABaseInteraction.h"
#include "DNA2ModInteraction.h"

/**
 * @brief CUDA implementation of the oxDNA model, as provided by DNAInteraction.
 */
template<typename number, typename number4>
class CUDADNA2ModInteraction: public CUDABaseInteraction<number, number4>, public DNA2ModInteraction<number> {
protected:
	// vectors that store the modification parameters for each nucleotide
	number *_d_stacking_roll, *_d_stacking_r_roll, *_d_stacking_tilt, *_d_stacking_multiplier, *_d_hb_multiplier;

public:
        enum {
		DEBYE_HUCKEL = 7
	};
	CUDADNA2ModInteraction();
	virtual ~CUDADNA2ModInteraction();

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

	void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_qorientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box);
};

extern "C" BaseInteraction<float> *make_CUDADNA2ModInteraction_float() { return new CUDADNA2ModInteraction<float, float4>(); } 
extern "C" BaseInteraction<double> *make_CUDADNA2ModInteraction_double() { return new CUDADNA2ModInteraction<double, LR_double4>(); }
#endif /* CUDADNA2MODINTERACTION_H_ */
