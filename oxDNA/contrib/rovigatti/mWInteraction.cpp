/*
 * mWInteraction.cpp
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#include "mWInteraction.h"
#include "../Particles/PatchyParticle.h"
#include "../Utilities/Utils.h"

#include <string>

using namespace std;

template <typename number>
mWInteraction<number>::mWInteraction() : BaseInteraction<number, mWInteraction<number> >(), _N(-1) {
	this->_int_map[mW] = &mWInteraction<number>::_two_body;

	_lambda = 1.;
	_gamma = 1.2;
	_a = 1.8;
	_A = 7.049556277;
	_B = 0.6022245584;
	_theta0 = 1.9106119321581925;
}

template <typename number>
mWInteraction<number>::~mWInteraction() {

}

template<typename number>
void mWInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	getInputNumber(&inp, "mW_lambda", &_lambda, 0);
	getInputNumber(&inp, "mW_gamma", &_gamma, 0);
	getInputNumber(&inp, "mW_theta0", &_theta0, 0);
	getInputNumber(&inp, "mW_a", &_a, 0);
	getInputNumber(&inp, "mW_A", &_A, 0);
	getInputNumber(&inp, "mW_B", &_B, 0);
}

template<typename number>
void mWInteraction<number>::init() {
	this->_rcut = _a;
	this->_sqr_rcut = SQR(this->_rcut);
	_cos_theta0 = cos(_theta0);

	OX_LOG(Logger::LOG_INFO, "mW parameters: lambda = %lf, A = %lf, B = %lf, gamma = %lf, theta0 = %lf, a = %lf", _lambda, _A, _B, _gamma, _theta0, _a);
}

template<typename number>
void mWInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) {
		particles[i] = new BaseParticle<number>();
	}
	_bonds.resize(N);
	_N = N;
}

template<typename number>
number mWInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number energy = pair_interaction_bonded(p, q, r, update_forces);
	energy += pair_interaction_nonbonded(p, q, r, update_forces);
	return energy;
}

template<typename number>
number mWInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	// reset _bonds. This is beyond horrible
	if(q == P_VIRTUAL && p->index == 0) {
		for(int i = 0; i < _N; i++)	_bonds[i].clear();
	}

	return (number) 0.f;
}

template<typename number>
number mWInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = q->pos.minimum_image(p->pos, this->_box_side);
		r = &computed_r;
	}

	return _two_body(p, q, r, update_forces);
}

template<typename number>
void mWInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	*N_strands = N;

	allocate_particles(particles, N);
	for (int i = 0; i < N; i ++) {
	   particles[i]->index = i;
	   particles[i]->type = P_A;
	   particles[i]->btype = P_A;
	   particles[i]->strand_id = i;
	}
}

template<typename number>
void mWInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

extern "C" mWInteraction<float> *make_float() {
	return new mWInteraction<float>();
}

extern "C" mWInteraction<double> *make_double() {
	return new mWInteraction<double>();
}

template class mWInteraction<float>;
template class mWInteraction<double>;
