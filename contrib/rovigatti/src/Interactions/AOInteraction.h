/*
 * AOInteraction.h
 *
 *  Created on: 24/oct/2017
 *      Author: lorenzo
 */

#ifndef AOINTERACTION_H_
#define AOINTERACTION_H_

#include "Interactions/BaseInteraction.h"

template <typename number>
class AOInteraction: public BaseInteraction<number, AOInteraction<number> > {
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

	virtual void allocate_particles(BaseParticle<number> **particles, int N);
	virtual void read_topology(int N, int *N_strands, BaseParticle<number> **particles);

	virtual number pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_term(int name, BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, r, update_forces);
	}

	virtual void check_input_sanity(BaseParticle<number> **particles, int N);
};

extern "C" AOInteraction<float> *make_AOInteraction_float();
extern "C" AOInteraction<double> *make_AOInteraction_double();

#endif /* AOINTERACTION_H_ */
