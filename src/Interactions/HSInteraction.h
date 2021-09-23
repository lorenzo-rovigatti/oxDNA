/*
 * HSInteraction.h
 *
 *  Created on: 28/Oct/2013
 *      Author: Flavio
 */

#ifndef HSINTERACTION_H_
#define HSINTERACTION_H_

#include "BaseInteraction.h"

/**
 * @brief Interaction class to simulate Hard spheres.
 *
 * Not much to discuss here. The radius is implicitely assued to bo 0.5
 * This interaction is selected with
 * interaction_type = HS
 */

class HSInteraction: public BaseInteraction {
protected:
	inline number _hs_pot(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

public:
	enum {
		HS = 0
	};

	HSInteraction();
	virtual ~HSInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
};

number HSInteraction::_hs_pot(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(update_forces) throw oxDNAException("No forces, figlio di ndrocchia");

	number rnorm = _computed_r.norm();
	number energy = 0;

	if(rnorm < (number) 1.) {
		this->set_is_infinite(true);
		//this->is_infinite = true;
		energy = (number) 1.e12;
	}
	else energy = (number) 0.;

	return energy;
}

#endif /* HSINTERACTION_H_ */

