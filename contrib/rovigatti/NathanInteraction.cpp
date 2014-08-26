/*
 * NathanInteraction.cpp
 *
 *  Created on: 21 Aug 2014
 *      Author: lorenzo
 */

#include "NathanInteraction.h"

template<typename number>
NathanInteraction<number>::NathanInteraction() {
	this->_int_map[PATCHY] = &NathanInteraction<number>::_patchy_interaction;

	_rep_E_cut = 0.;
	_rep_power = 200;

	_patch_alpha = 0.12;
	_patch_power = 30;
}

template<typename number>
NathanInteraction<number>::~NathanInteraction() {

}

template<typename number>
void NathanInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	BaseInteraction<number, NathanInteraction<number> >::read_topology(N, N_strands, particles);
}

template<typename number>
void NathanInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) {
		particles[i] = new NathanPatchyParticle<number>();
	}
}

template<typename number>
void NathanInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);

	getInputNumber(&inp, "NATHAN_alpha", &_patch_alpha, 0);
	getInputNumber(&inp, "NATHAN_cosmax", &_patch_cosmax, 1);
}

template<typename number>
void NathanInteraction<number>::init() {
	_patch_cutoff = _patch_alpha * 1.5;
	this->_rcut = 1. + _patch_cutoff;
	this->_sqr_rcut = SQR(this->_rcut);

	_rep_E_cut = pow((number) this->_rcut, -_rep_power);

	_patch_pow_sigma = pow(_patch_cosmax, _patch_power);
	_patch_pow_alpha = pow(_patch_alpha, (number) 10.);
	number r8b10 = pow(_patch_cutoff, (number) 8.) / _patch_pow_alpha;
	_patch_E_cut = -1.001 * exp(-(number)0.5 * r8b10 * SQR(_patch_cutoff));
	_patch_E_cut = 0.;
}

template<typename number>
void NathanInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template<typename number>
number NathanInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number NathanInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return (number) 0.f;
}

template<typename number>
number NathanInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = q->pos.minimum_image(p->pos, this->_box_side);
		r = &computed_r;
	}

	return _patchy_interaction(p, q, r, update_forces);
}

template class NathanInteraction<float>;
template class NathanInteraction<double>;
