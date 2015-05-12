/*
 * GraftedInteraction.h
 *
 *  Created on: 10/feb/2013
 *      Author: lorenzo
 */

#ifndef GRAFTEDINTERACTION_H_
#define GRAFTEDINTERACTION_H_

#define P_COLLOID 2

#include "BaseInteraction.h"
#include "TSPInteraction.h"

/**
 * @brief Handles the interaction between coarse-grained DNA tetramers.
 */
template <typename number>
class GraftedInteraction: public BaseInteraction<number, GraftedInteraction<number> > {
protected:
	TSPInteraction<number> _TSP_inter;

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

	virtual void set_box_side(number box_side);
	virtual void set_box(BaseBox<number> *box);

	virtual void allocate_particles(BaseParticle<number> **particles, int N);
	virtual void read_topology(int N, int *N_strands, BaseParticle<number> **particles);

	virtual number _wall_interaction(BaseParticle<number> *p, bool update_forces);

	virtual number pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_term(int name, BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, r, update_forces);
	}

	virtual void check_input_sanity(BaseParticle<number> **particles, int N);
	virtual void generate_random_configuration(BaseParticle<number> **particles, int N, number box_side);
};

template<typename number>
class Colloid: public TSPParticle<number> {
protected:
	std::vector<LR_vector<number> > _base_anchors;

public:
	Colloid(std::vector<LR_vector<number> > &anchor_poss);
	virtual ~Colloid();

	virtual bool is_rigid_body() { return true; }

	void set_positions();
};

extern "C" GraftedInteraction<float> *make_float() { return new GraftedInteraction<float>(); }
extern "C" GraftedInteraction<double> *make_double() { return new GraftedInteraction<double>(); }

#endif /* GRAFTEDINTERACTION_H_ */
