/*
 * AOInteraction.h
 *
 *  Created on: 24/oct/2017
 *      Author: lorenzo
 */

#ifndef AOINTERACTION_H_
#define AOINTERACTION_H_

#include "Interactions/BaseInteraction.h"

class AOInteraction: public BaseInteraction {
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

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);
	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
};

extern "C" AOInteraction *make_AOInteraction();

#endif /* AOINTERACTION_H_ */
