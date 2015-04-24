
#ifndef TEP_INTERACTION_H
#define TEP_INTERACTION_H

#include "BaseInteraction.h"
#include "../Particles/TEPParticle.h"

/**
 * @brief Handles interactions between elements of the Twistable Elastic Polymer (TEP) model.
 *
 * The model is described in the paper (add reference), eq.4.
 * This description should be edited when the model has finally been developed.
 *
 * Input options:
 *
 * @verbatim
[_prefer_harmonic_over_fene = <bool>(if True, neighbouring beads are bound by an harmonic potential instead of a FENE one. Defaults to false.)]
[_allow_broken_fene = <bool>(if True, the code won't die when the beads are beyond the acceptable FENE range. Defaults to True.)]
@endverbatim
 The option above are true for the DNAinteraction class, not this one, but I keep them so that
 I can keep track of how they should be declared.
 */
template <typename number>
class TEPInteraction : public BaseInteraction<number, TEPInteraction<number> > {
protected:
	number _T;
	//parameters of the TEP model
	number _ka;
	number _kb;
	number _kt;
	// FENE parameters
	number _TEP_FENE_DELTA;
	number _TEP_FENE_DELTA2;
	number _TEP_FENE_EPS;
	number _TEP_FENE_R0;
	// LJ parameters
	number _TEP_EXCL_EPS_BONDED;
	number _TEP_EXCL_EPS_NONBONDED;
	number _TEP_EXCL_S2;
	number _TEP_EXCL_R2;
	number _TEP_EXCL_B2;
	number _TEP_EXCL_RC2;
	
	number _TEP_spring_offset;
	/// true by default; set this to false if you want the code to not die when bonded backbones are found to be outside the acceptable FENE range
	bool _allow_broken_fene;

// false by default: set this to true if you want neighbouring particles to be bound by a quadratic potential instead of a FENE
	bool _prefer_harmonic_over_fene;



	virtual number _spring(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _bonded_bending(BaseParticle <number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _bonded_twist(BaseParticle <number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _bonded_alignment(BaseParticle <number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _bonded_debye_huckel(BaseParticle <number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _nonbonded_debye_huckel(BaseParticle <number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _bonded_excluded_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _nonbonded_excluded_volume(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);

	inline number _repulsive_lj2(number prefactor,const LR_vector<number> &r, LR_vector<number> &force, number sigma, number rstar, number b, number rc, bool update_forces);
	void setNonNegativeNumber(input_file *inp, const char * skey, number *dest, int mandatory, const char * arg_description);
	void setPositiveNumber(input_file *inp, const char * skey, number *dest, int mandatory, const char * arg_description);


	/**
	 * @brief Checks whether the two particles share a backbone link.
	 *
	 * @param p
	 * @param q
	 * @return true if they are bonded, false otherwise
	 */
	inline bool _are_bonded(BaseParticle<number> *p, BaseParticle<number> *q) { return (p->n3 == q || p->n5 == q); }

public:
	enum {
		SPRING = 0,
		BONDED_BENDING = 1,
		BONDED_TWIST = 2,
		BONDED_ALIGNMENT = 3,
		NONBONDED_EXCLUDED_VOLUME = 4,
		BONDED_DEBYE_HUCKEL = 5,
		NONBONDED_DEBYE_HUCKEL = 6
	};
	TEPInteraction();
	virtual ~TEPInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(BaseParticle<number> **particles, int N);

	virtual number pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_term(int name, BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, r, update_forces);
	}

	virtual void check_input_sanity(BaseParticle<number> **particles, int N);

	virtual void read_topology(int N_from_conf, int *N_strands, BaseParticle<number> **particles);

};

template<typename number>
number TEPInteraction<number>::_repulsive_lj2(number prefactor, const LR_vector<number> &r, LR_vector<number> &force, number sigma, number rstar, number b, number rc, bool update_forces) {
	// this is a bit faster than calling r.norm()
	number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
	number energy = (number) 0;
	if(rnorm < SQR(rc)) {
		if(rnorm > SQR(rstar)) {
			number rmod = sqrt(rnorm);
			number rrc = rmod - rc;
			energy = prefactor * b * SQR(rrc);
			if(update_forces) force = -r * (2 * prefactor * b * rrc / rmod);
		}
		else {
			number tmp = SQR(sigma) / rnorm;
			number lj_part = tmp * tmp * tmp;
			// the additive term was added by me to mimick Davide's implementation
			energy = 4 * prefactor * (SQR(lj_part) - lj_part) + prefactor; 
			if(update_forces) force = -r * (24 * prefactor * (lj_part - 2*SQR(lj_part)) / rnorm);
		}
	}

	if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;

	return energy;
}

template<typename number>
	void TEPInteraction<number>::setNonNegativeNumber(input_file *inp, const char * skey, number *dest, int mandatory, const char * arg_description){

	if(getInputNumber(inp, skey,dest, mandatory) == KEY_FOUND) {
		if( *dest < 0 ) throw oxDNAException("read negative parameter %s (%s) for the TEP model. %s = %g. Aborting",skey,arg_description,skey,*dest); 
		OX_LOG(Logger::LOG_INFO,"%s manually set to %g",skey,*dest);
	}
}

template<typename number>
	void TEPInteraction<number>::setPositiveNumber(input_file *inp, const char * skey, number *dest, int mandatory, const char * arg_description){

	if(getInputNumber(inp, skey,dest, mandatory) == KEY_FOUND) {
		if( *dest <= 0 ) throw oxDNAException("read non-positive parameter %s (%s) for the TEP model. %s = %g. Aborting",skey,arg_description,skey,*dest); 
		OX_LOG(Logger::LOG_INFO,"%s manually set to %g",skey,*dest);
	}
}

#endif /* TEP_INTERACTION_H */
