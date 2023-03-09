/**
 * @brief Enables simulation of DNA-RNA hybrid systems with the DNA2 and RNA2 interaction types.
 * 
 * To use the model, set interaction_type = DNA2withRNA2 in the input file
 *
 * 
 * 
 * NB: 
 * The debye-huckel interaction used is copied from DNA2, and is the same for all interaction
 * types in the system (DNA-DNA, RNA-RNA and hybrid). The prefactor has been set to 0.07, which is 
 * the mean of the default values used for DNA and RNA alone.
 * 		
 * Other things which haven't been included:
 * -RNA mismatch potential
 * -Coaxial stacking between DNA and RNA
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * Input options:
 * 
 * @verbatim
 salt_concentration = <float>  (sets the salt concentration in M)
 [dh_lambda = <float> (the value that lambda, which is a function of temperature (T) and salt concentration (I), should take when T=300K and I=1M, defaults to the value from Debye-Huckel theory, 0.3616455)]
 [dh_strength = <float> (the value that scales the overall strength of the Debye-Huckel interaction, defaults to 0.0543)]
 [dh_half_charged_ends = <bool>  (set to false for 2N charges for an N-base-pair duplex, defaults to 1)]
 @endverbatim
 */

#ifndef DNA2_WITH_RNA2_INTERACTION_H
#define DNA2_WITH_RNA2_INTERACTION_H

#include "DNAwithRNAInteraction.h"

class DNA2withRNA2Interaction: public DNAwithRNAInteraction {

protected:
	float _salt_concentration;
	bool _debye_huckel_half_charged_ends;
	number _debye_huckel_prefactor;
	number _debye_huckel_lambdafactor;

	//the following values are calculated
	number _debye_huckel_RC; // this is the maximum interaction distance between backbones to interact with DH
	number _debye_huckel_RHIGH; // distance after which the potential is replaced by a quadratic cut-off
	number _debye_huckel_B; // prefactor of the quadratic cut-off
	number _minus_kappa;

	number _f4_pure_harmonic_DNA(number t, int type);
	number _f4Dsin_pure_harmonic_DNA(number t, int type);
	number _f4D_pure_harmonic_DNA(number t, int type);

	/**
	 * @brief Custom function that returns f4. Meshes are never used for the coaxial stacking f4(theta1) term for the DNA2 interaction
	 *
	 * @param cost  argument of f4
	 * @param i     type of the interaction (which mesh to use)
	 */
	virtual number _custom_f4_DNA(number cost, int i) {
		if(i != CXST_F4_THETA1) return this->_mesh_f4_DNA[i].query(cost);
		else return this->_fakef4_cxst_t1_DNA(cost, (void *) &i);
	}

	/**
	 * @brief Custom function that returns the derivative of f4. See _custom_f4_DNA
	 *
	 * @param cost  argument of f4D
	 * @param i     type of the interaction (which mesh to use)
	 */
	virtual number _custom_f4D_DNA(number cost, int i) {
		if(i != CXST_F4_THETA1) return this->_mesh_f4_DNA[i].query_derivative(cost);
		else return this->_fakef4D_cxst_t1_DNA(cost, (void *) &i);
	}

public:
	virtual number _debye_huckel(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	enum {
		DEBYE_HUCKEL = 7
	};
	DNA2withRNA2Interaction();
	virtual ~DNA2withRNA2Interaction() {
	}

	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	virtual void get_settings(input_file &inp);
	virtual void init();

	number _fakef4_cxst_t1_DNA(number t, void * par);
	number _fakef4D_cxst_t1_DNA(number t, void * par);

	virtual number _coaxial_stacking_DNA(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	number F4_THETA_SA_DNA[13];
	number F4_THETA_SB_DNA[13];
};




#endif /* DNA2_INTERACTION_H */
