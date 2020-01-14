/*
 * mWInteraction.h
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#ifndef MWINTERACTION_H_
#define MWINTERACTION_H_

#include "Interactions/BaseInteraction.h"

struct mWBond {
	BaseParticle *other;
	LR_vector r;
	number mod_r;

	mWBond(BaseParticle *o, LR_vector my_r, number my_mod_r) :
					other(o),
					r(my_r),
					mod_r(my_mod_r) {
	}
};

/**
 * @brief Manages the interaction between mW-like patchy particles.
 *
 * This interaction is selected with
 * interaction_type = patchy
 *
 * @verbatim
 mW_N = <int> (number of patches)
 [mW_N_B = <int> (number of patches on species B)]
 [mW_one_component = <bool> (if true then the system will not contain particles of type B and A-A bonds will be allowed. Defaults to false.)]
 @endverbatim
 */
class mWInteraction: public BaseInteraction<mWInteraction> {
protected:
	int _N;

	number _A, _B;
	number _a;
	number _gamma;
	number _lambda;
	number _theta0, _cos_theta0;

	std::vector<std::vector<mWBond> > _bonds;

	/**
	 * @brief mW interaction between two particles.
	 *
	 * @param p
	 * @param q
	 * @param r
	 * @param update_forces
	 * @return
	 */
	inline number _two_body(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces);

	inline number _three_body(BaseParticle *p, mWBond &new_bond, bool update_forces);

public:
	enum {
		mW = 0
	};

	mWInteraction();
	virtual ~mWInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles, int N);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false);
	virtual number pair_interaction_term(int name, BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, r, update_forces);
	}

	virtual void read_topology(int N, int *N_strands, std::vector<BaseParticle *> &particles);
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles, int N);
};

number mWInteraction::_two_body(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {
	number sqr_r = r->norm();
	if(sqr_r > this->_sqr_rcut) return (number) 0.f;

	number energy = (number) 0.f;

	// centre-centre
	if(sqr_r < this->_sqr_rcut) {
		number ir4 = 1. / SQR(sqr_r);
		number mod_r = sqrt(sqr_r);
		number mod_r_a = mod_r - _a;
		number exp_part = exp(1. / mod_r_a);
		energy = _A * (_B * ir4 - 1.) * exp_part;

		mWBond p_bond(q, *r, mod_r);
		mWBond q_bond(p, -(*r), mod_r);

		if(update_forces) {
			LR_vector force = *r * ((_A * 4. * exp_part * _B * ir4 / mod_r + energy / SQR(mod_r_a)) / mod_r);
			p->force -= force;
			q->force += force;
		}

		energy += _three_body(p, p_bond, update_forces);
		energy += _three_body(q, q_bond, update_forces);

		_bonds[p->index].push_back(p_bond);
		_bonds[q->index].push_back(q_bond);
	}

	return energy;
}

number mWInteraction::_three_body(BaseParticle *p, mWBond &new_bond, bool update_forces) {
	number energy = 0.;

	typename std::vector<mWBond>::iterator it = _bonds[p->index].begin();
	for(; it != _bonds[p->index].end(); it++) {
		number irpq = it->mod_r * new_bond.mod_r;
		number cos_theta = it->r * new_bond.r / irpq;
		number diff_cos = cos_theta - _cos_theta0;
		number exp_part = exp(_gamma / (it->mod_r - _a)) * exp(_gamma / (new_bond.mod_r - _a));
		number l_diff_exp = _lambda * diff_cos * exp_part;
		number U3 = l_diff_exp * diff_cos;
		energy += U3;

		if(update_forces) {
			LR_vector p_it_force = it->r * (U3 * _gamma / SQR(it->mod_r - _a) / it->mod_r + 2. * cos_theta * l_diff_exp / SQR(it->mod_r)) - new_bond.r * (2. * l_diff_exp / irpq);
			p->force -= p_it_force;
			it->other->force += p_it_force;

			LR_vector q_it_force = new_bond.r * (U3 * _gamma / SQR(new_bond.mod_r - _a) / new_bond.mod_r + 2. * cos_theta * l_diff_exp / SQR(new_bond.mod_r)) - it->r * (2. * l_diff_exp / irpq);
			p->force -= q_it_force;
			new_bond.other->force += q_it_force;
		}
	}

	return energy;
}

extern "C" mWInteraction *make_mWInteraction();

#endif /* MWINTERACTION_H_ */
