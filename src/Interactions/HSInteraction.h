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
template <typename number>
class HSInteraction: public BaseInteraction<number, HSInteraction<number> > {
protected:
	inline number _hs_pot (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);

public:
	enum {
		HS = 0
	};

	HSInteraction();
	virtual ~HSInteraction();

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

	//virtual void generate_random_configuration(BaseParticle<number> **particles, int N, number box_side);
};

template<typename number>
number HSInteraction<number>::_hs_pot (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if (update_forces) throw oxDNAException ("No forces, figlio di ndrocchia");

	number rnorm = r->norm();
	number energy = 0;

	if (rnorm < (number) 1.) {
		this->set_is_infinite (true);
		//this->is_infinite = true;
		energy = (number) 1.e12;
	}
	else energy = (number) 0.;

	return energy;
}


#endif /* HSINTERACTION_H_ */

