/*
 * GraftedInteraction.h
 *
 *  Created on: 10/feb/2013
 *      Author: lorenzo
 */

#ifndef GRAFTEDINTERACTION_H_
#define GRAFTEDINTERACTION_H_

#define P_COLLOID 2

#include "Interactions/BaseInteraction.h"
#include "TSPInteraction.h"

/**
 * @brief Handles the interaction between coarse-grained DNA tetramers.
 */
class GraftedInteraction: public BaseInteraction {
protected:
	TSPInteraction _TSP_inter;

	int _N_arms, _N_per_arm;
	number _alpha;
	number _colloid_sigma;
	number _colloid_monomer_sqr_sigma;
	number _colloid_monomer_sqr_rep_rcut;
	number _colloid_rfene, _colloid_sqr_rfene;
	int _colloid_n;

	bool _walls;
	number _wall_distance;

	char _anchor_file[512];

public:
	GraftedInteraction();
	virtual ~GraftedInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void set_box(BaseBox *box);

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);
	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);

	virtual number _wall_interaction(BaseParticle *p, bool update_forces);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
	virtual void generate_random_configuration(std::vector<BaseParticle *> &particles);
};

class Colloid: public TSPParticle {
protected:
	std::vector<LR_vector> _base_anchors;

public:
	Colloid(std::vector<LR_vector> &anchor_poss);
	virtual ~Colloid();

	virtual bool is_rigid_body() {
		return true;
	}

	void set_positions();
};

extern "C" GraftedInteraction *make_GraftedInteraction() {
	return new GraftedInteraction();
}

#endif /* GRAFTEDINTERACTION_H_ */
