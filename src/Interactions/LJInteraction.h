/*
 * LJInteraction.h
 *
 *  Created on: 10/feb/2013
 *      Author: lorenzo
 */

#ifndef LJINTERACTION_H_
#define LJINTERACTION_H_

#include "BaseInteraction.h"

/**
 * @brief Handles (generalised) Lennard-Jones interactions between spheres of size 1 or a Kob-Andersen interaction.
 *
 * TODO: the Kob-Andersen mixture should be better implemented.
 *
 * This interaction is selected with
 * interaction_type = LJ
 *
 * Input options:
 *
 * @verbatim
LJ_rcut = <float> (interaction cutoff)
[LJ_kob_andersen = <bool> (Simulate a Kob-Andersen mixture. Defaults to false.)]
[LJ_n = <int> (Generalised LJ exponent. Defaults to 6, which is the classic LJ value.)]
@endverbatim
 */
template <typename number>
class LJInteraction: public BaseInteraction<number, LJInteraction<number> > {
protected:
	number _E_cut[3];
	bool _is_ka_mixture;
	int _N_A, _N_B;
	int _n[3];
	number _sigma[3];
	number _sqr_sigma[3];
	number _epsilon[3];
	number _sqr_LJ_rcut[3];

	inline number _lennard_jones(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);

public:
	enum {
		LENNARD_JONES = 0
	};

	LJInteraction();
	virtual ~LJInteraction();

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

template<typename number>
number LJInteraction<number>::_lennard_jones(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number rnorm = r->norm();
	number energy = 0;

	int type = p->type + q->type;
	if(rnorm < _sqr_LJ_rcut[type]) {
		number tmp = this->_sqr_sigma[type] / rnorm;
		number lj_part = 1;
		for(int i = 0; i < _n[type]/2; i++) lj_part *= tmp;
		energy = 4 * _epsilon[type] * (SQR(lj_part) - lj_part) - _E_cut[type];
		if(update_forces) {
			LR_vector<number> force = *r * (-4 *_n[type] * _epsilon[type] * (lj_part - 2*SQR(lj_part)) / rnorm);
			p->force -= force;
			q->force += force;
		}
	}

	return energy;
}

#endif /* LJINTERACTION_H_ */
