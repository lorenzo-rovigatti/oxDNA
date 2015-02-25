/**
 * @brief Adds a Debye-Huckel (salt dependent) term to the oxDNA potential represented in DNAInteraction
 *
 * Ben Jan 2015
 * Ferdinando (April 2014) after petr (April 2014)
 * 
 * To use the oxDNA2 model, set
 *
 * interaction_type = DNA2 
 *
 * in the input file
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

#ifndef DNA2_INTERACTION_H
#define DNA2_INTERACTION_H

#include "DNAInteraction.h"

template<typename number>
class DNA2Interaction: public DNAInteraction<number> {

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

	number _f4_pure_harmonic(number t, int type);
	number _f4Dsin_pure_harmonic(number t, int type);
	number _f4D_pure_harmonic(number t, int type);

	/**
	 * @brief Custom function that returns f4. Meshes are never used for the coaxial stacking f4(theta1) term for the DNA2 interaction
	 *
	 * @param cost  argument of f4
	 * @param i     type of the interaction (which mesh to use)
	 */
	virtual number _custom_f4 (number cost, int i) { 
		if (i != CXST_F4_THETA1) return this->_query_mesh (cost, this->_mesh_f4[i]); 
		else return this->_fakef4_cxst_t1 (cost, (void *)&i);
	}

	/**
	 * @brief Custom function that returns the derivative of f4. See _custom_f4
	 *
	 * @param cost  argument of f4D
	 * @param i     type of the interaction (which mesh to use)
	 */
	virtual number _custom_f4D (number cost, int i) { 
		if (i != CXST_F4_THETA1) return this->_query_meshD (cost, this->_mesh_f4[i]); 
		else return this->_fakef4D_cxst_t1 (cost, (void *)&i);
	}

public:
	virtual number _debye_huckel(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);

        enum {
		DEBYE_HUCKEL = 7
	};
	DNA2Interaction();
	virtual ~DNA2Interaction() {}
    
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);

	virtual void check_input_sanity(BaseParticle<number> **particles, int N);

	virtual void get_settings(input_file &inp);
	virtual void init();

	number _fakef4_cxst_t1(number t, void * par);
	number _fakef4D_cxst_t1(number t, void * par);

	virtual number _backbone(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _coaxial_stacking(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);

	number F4_THETA_SA[13];
	number F4_THETA_SB[13];
};

/**
 * @brief Handles interactions between DNA nucleotides without using meshes.
 *
 * This interaction is selected with
 * interaction_type = DNA2_nomesh
 */
template <typename number>
class DNA2Interaction_nomesh : public DNA2Interaction<number> {
protected:

	/**
	 * @brief Custom function that returns f4.
	 *
	 * @param cost  argument of f4
	 * @param i     type of the interaction
	 */
	virtual number _custom_f4 (number cost, int i) {
		if (i != CXST_F4_THETA1) return this->_fakef4 (cost, (void *)&i);
		else return this->_fakef4_cxst_t1 (cost, (void *)&i);
	}

	/**
	 * @brief Custom function that returns the derivative of f4. See _custom_f4
	 *
	 * @param cost  argument of f4D
	 * @param i     type of the interaction
	 */
	virtual number _custom_f4D (number cost, int i) {
		if (i != CXST_F4_THETA1) return this->_fakef4D (cost, (void *)&i);
		else return this->_fakef4D_cxst_t1 (cost, (void *)&i);
	}

public:
	DNA2Interaction_nomesh(){};
	virtual ~DNA2Interaction_nomesh(){};

};

template class DNA2Interaction_nomesh<float>;
template class DNA2Interaction_nomesh<double>;

#endif /* DNA2_INTERACTION_H */
