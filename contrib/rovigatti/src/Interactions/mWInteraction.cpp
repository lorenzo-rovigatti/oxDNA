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
				BaseInteraction<mWInteraction>(),
				_N(-1) {
	_int_map[mW] = &mWInteraction::_two_body;

	_lambda = 1.;
	_gamma = 1.2;
	_a = 1.8;
	_A = 7.049556277;
	_B = 0.6022245584;
	_theta0 = 1.9106119321581925;
}

mWInteraction::~mWInteraction() {

}

void mWInteraction::get_settings(input_file &inp) {
	IBaseInteraction::get_settings(inp);

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
