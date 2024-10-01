/*
 * mWInteraction.cpp
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#include "mWInteraction.h"
#include "Utilities/Utils.h"

#include <string>

using namespace std;

mWInteraction::mWInteraction() :
				BaseInteraction(),
				_N(-1) {
	ADD_INTERACTION_TO_MAP(mW, _two_body);

	_lambda = 1.;
	_gamma = 1.2;
	_a = 1.8;
	_A = 7.049556277;
	_B = 0.6022245584;
	_theta0 = 1.9106119321581925;
	_cos_theta0 = cos(_theta0);
}

mWInteraction::~mWInteraction() {

}

void mWInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	getInputNumber(&inp, "mW_lambda", &_lambda, 0);
	getInputNumber(&inp, "mW_gamma", &_gamma, 0);
	getInputNumber(&inp, "mW_theta0", &_theta0, 0);
	getInputNumber(&inp, "mW_a", &_a, 0);
	getInputNumber(&inp, "mW_A", &_A, 0);
	getInputNumber(&inp, "mW_B", &_B, 0);
}

void mWInteraction::init() {
	_rcut = _a;
	_sqr_rcut = SQR(_rcut);
	_cos_theta0 = cos(_theta0);

	OX_LOG(Logger::LOG_INFO, "mW parameters: lambda = %lf, A = %lf, B = %lf, gamma = %lf, theta0 = %lf, a = %lf", _lambda, _A, _B, _gamma, _theta0, _a);
}

void mWInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
	int N = particles.size();
	for(int i = 0; i < N; i++) {
		particles[i] = new BaseParticle();
	}
	_bonds.resize(N);
	_N = N;
}

void mWInteraction::begin_energy_computation() {
	BaseInteraction::begin_energy_computation();

	for(int i = 0; i < _N; i++) {
		_bonds[i].clear();
	}
}

number mWInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		if(q != P_VIRTUAL && p != P_VIRTUAL) {
			_computed_r = _box->min_image(p->pos, q->pos);
		}
	}

	number energy = pair_interaction_bonded(p, q, false, update_forces);
	energy += pair_interaction_nonbonded(p, q, false, update_forces);
	return energy;
}

number mWInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return (number) 0.f;
}

number mWInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	return _two_body(p, q, false, update_forces);
}

number mWInteraction::_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm();
	if(sqr_r > _sqr_rcut) {
		return (number) 0.f;
	}

	number energy = (number) 0.f;

	// centre-centre
	if(sqr_r < _sqr_rcut) {
		number ir4 = 1. / SQR(sqr_r);
		number mod_r = sqrt(sqr_r);
		number mod_r_a = mod_r - _a;
		number exp_part = exp(1. / mod_r_a);
		energy = _A * (_B * ir4 - 1.) * exp_part;

		mWBond p_bond(q, _computed_r, mod_r);
		mWBond q_bond(p, -_computed_r, mod_r);

		if(update_forces) {
			LR_vector force = _computed_r * ((_A * 4. * exp_part * _B * ir4 / mod_r + energy / SQR(mod_r_a)) / mod_r);
			p->force -= force;
			q->force += force;

			_update_stress_tensor(p->pos, -force);
			_update_stress_tensor(p->pos + _computed_r, force);
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

			_update_stress_tensor(p->pos + it->r, p_it_force);
			_update_stress_tensor(p->pos + new_bond.r, q_it_force);
			_update_stress_tensor(p->pos, -(p_it_force + q_it_force));
		}
	}

	return energy;
}

void mWInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
	int N = particles.size();
	*N_strands = N;

	allocate_particles(particles);
	for(int i = 0; i < N; i++) {
		particles[i]->index = i;
		particles[i]->type = P_A;
		particles[i]->btype = P_A;
		particles[i]->strand_id = i;
	}
}

void mWInteraction::check_input_sanity(std::vector<BaseParticle *> &particles) {

}

extern "C" mWInteraction *make_mWInteraction() {
	return new mWInteraction();
}
