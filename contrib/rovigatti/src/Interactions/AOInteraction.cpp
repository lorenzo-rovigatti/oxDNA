/*
 * AOInteraction.cpp
 *
 *  Created on: 24/oct/2017
 *      Author: lorenzo
 */

#include "AOInteraction.h"

AOInteraction::AOInteraction() :
				BaseInteraction() {
	ADD_INTERACTION_TO_MAP(AO, pair_interaction_nonbonded);

	_colloid_sigma_sqr = 1. / SQR(1.01557);
	_h_zausch_4 = SQR(SQR(0.01)*_colloid_sigma_sqr);
	_rep_rcut_sqr = pow(2., 1. / 3.) * _colloid_sigma_sqr;
	_rep_rcut = sqrt(_rep_rcut_sqr);
}

AOInteraction::~AOInteraction() {

}

void AOInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	getInputNumber(&inp, "AO_strength", &_attraction_strength, 1);
	getInputNumber(&inp, "AO_q", &_q, 1);
}

void AOInteraction::init() {
	_sigma_colloid_polymer = (1 + _q) / 2.;
	this->_rcut = 1 + _q;
	if(this->_rcut < _rep_rcut) {
		throw oxDNAException("AO: 1 + q (%lf) is smaller than rep_rcut (%lf): this is not supported", 1 + _q, _rep_rcut);
	}
	this->_sqr_rcut = SQR(this->_rcut);

	OX_LOG(Logger::LOG_INFO, "AO: rcut = %lf, rep_rcut = %lf", this->_rcut, _rep_rcut);
}

void AOInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
	for(uint i = 0; i < particles.size(); i++) {
		particles[i] = new BaseParticle();
	}
}

void AOInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
	int N = particles.size();
	*N_strands = N;

	allocate_particles(particles);
	for(int i = 0; i < N; i++) {
		particles[i]->index = particles[i]->strand_id = i;
		particles[i]->type = particles[i]->btype = P_A;
	}
}

number AOInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, compute_r, update_forces);
}

number AOInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return (number) 0.f;
}

number AOInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = this->_box->min_image(p->pos, q->pos);
	}

	number r_norm = _computed_r.norm();
	number energy = 0;

	if(r_norm < this->_sqr_rcut) {
		number r_mod = sqrt(r_norm);
		if(r_norm < _rep_rcut_sqr) {
			number WCA_part = CUB(_colloid_sigma_sqr / r_norm);
			number WCA = 4 * (SQR(WCA_part) - WCA_part + 0.25);
			number S_part = SQR(SQR(r_mod - _rep_rcut));
			number S = S_part / (_h_zausch_4 + S_part);
			energy += WCA * S;

			if(update_forces) {
				number WCA_der = 24. * (WCA_part - 2 * SQR(WCA_part)) / r_mod;
				number S_der = (1 - S) * (4 * CUB(r_mod - _rep_rcut)) / (_h_zausch_4 + S_part);
				number force = -(WCA_der * S + WCA * S_der);
				p->force -= _computed_r * (force / r_mod);
				q->force += _computed_r * (force / r_mod);
			}
		}

		number r_rescaled = r_mod / _sigma_colloid_polymer;
		energy += -_attraction_strength * (1. - 3. / 4. * r_rescaled + CUB(r_rescaled) / 16.);

		if(update_forces) {
			number force = -_attraction_strength * (3. / 4. - 3. * SQR(r_rescaled) / 16.) / _sigma_colloid_polymer;
			p->force -= _computed_r * (force / r_mod);
			q->force += _computed_r * (force / r_mod);
		}
	}

	return energy;
}

void AOInteraction::check_input_sanity(std::vector<BaseParticle *> &particles) {

}

extern "C" AOInteraction *make_AOInteraction() {
	return new AOInteraction();
}
