/*
 * AOInteraction.h
 *
 *  Created on: 24/oct/2017
 *      Author: lorenzo
 */

#ifndef AOINTERACTION_H_
#define AOINTERACTION_H_

#include "Interactions/BaseInteraction.h"

class AOInteraction: public BaseInteraction<AOInteraction> {
protected:
	number _q, _sigma_colloid_polymer;
	number _attraction_strength;
	number _colloid_sigma_sqr;
	number _h_zausch_4;
	number _rep_rcut, _rep_rcut_sqr;

public:
	enum {
		AO = 4
	};

	AOInteraction();
	virtual ~AOInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(BaseParticle **particles, int N);
	virtual void read_topology(int N, int *N_strands, BaseParticle **particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false);
	virtual number pair_interaction_term(int name, BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, r, update_forces);
	}

	virtual void check_input_sanity(BaseParticle **particles, int N);
};

extern "C" AOInteraction *make_AOInteraction();

#endif /* AOINTERACTION_H_ */
