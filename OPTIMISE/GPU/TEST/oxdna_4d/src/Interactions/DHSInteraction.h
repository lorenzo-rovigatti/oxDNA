/*
 * DHSInteraction.h
 *
 *  Created on: 28/Oct/2013
 *      Author: Flavio
 */

#ifndef DHSINTERACTION_H_
#define DHSINTERACTION_H_

#include "BaseInteraction.h"

/**
 * @brief Interaction class to simulate dipolar hard sphere
 * 
 * The dipolar interactions are simulated with a reaciot field. The reaction
 * field parameters can be input from the input file
 *
 * interaction_type = DHS
 * Input options:
 *
 * @verbatim
 DHS_eps = <float> (background dielectrci constant for reaction field treatment)
 DHS_rcut = <float> (cutoff for the reaction field treatment)
 @endverbatim
 */

class DHSInteraction: public BaseInteraction {
protected:
	number _eps;
	number _rf_fact;

	inline number _dhs_pot(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

public:
	enum {
		DHS = 0
	};

	DHSInteraction();
	virtual ~DHSInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);

	//virtual void generate_random_configuration(std::vector<BaseParticle *> &particles, number box_side);
};

number DHSInteraction::_dhs_pot(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(update_forces) {
		throw oxDNAException("No forces, figlio di ndrocchia");
	}

	number rnorm = _computed_r.norm();

	if(rnorm < (number) 1.) {
		this->set_is_infinite(true);
		return (number) 1.e12;
	}

	if(rnorm > this->_sqr_rcut) return (number) 0.f;

	// if things do not overlap, compute the DHS part of the
	// interaction potential
	number rnorm1 = sqrt(rnorm);
	number dot = p->orientation.v1 * q->orientation.v1;

	// direct part
	number energy = -(1. / (rnorm * rnorm1)) * (3 * (p->orientation.v1 * _computed_r) * (q->orientation.v1 * _computed_r) / rnorm - dot);

	// reaction field part
	energy += -_rf_fact * dot;

	return energy;
}

#endif /* HSINTERACTION_H_ */

