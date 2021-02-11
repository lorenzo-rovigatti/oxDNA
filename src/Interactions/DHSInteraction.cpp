/*
 * DHSInteraction.cpp
 *
 *  Created on: 29/Oct/2013
 *      Author: Flavio 
 */

#include "DHSInteraction.h"

DHSInteraction::DHSInteraction() :
				BaseInteraction() {
	ADD_INTERACTION_TO_MAP(DHS, _dhs_pot);
}

DHSInteraction::~DHSInteraction() {

}

void DHSInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);
	char tmps[512];
	getInputString(&inp, "sim_type", (char *) tmps, 1);
	if(strncmp(tmps, "MC", 512)) throw oxDNAException("Cannot run DHS with MD");

	// rcut for dipolar (and Reaction Field)
	float tmpf;
	getInputFloat(&inp, "DHS_rcut", &tmpf, 1);
	_rcut = (number) tmpf;

	// Reaction field medium epsilon
	getInputFloat(&inp, "DHS_eps", &tmpf, 1);
	_eps = (number) tmpf;

	_rf_fact = (_eps - 1.) / (2. * _eps + 1.) / (_rcut * _rcut * _rcut);

	OX_LOG(Logger::LOG_INFO, "Initializing Dipolar Hard Sphere potential with rcut = %g and dielectric constant %g", _rcut, _rf_fact);
}

void DHSInteraction::init() {
	_sqr_rcut = SQR(_rcut);
}

void DHSInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
	for(uint i = 0; i < particles.size(); i++) {
		particles[i] = new BaseParticle();
	}
}

number DHSInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, compute_r, update_forces);
}

number DHSInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return (number) 0.f;
}

number DHSInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	return _dhs_pot(p, q, false, update_forces);
}

void DHSInteraction::check_input_sanity(std::vector<BaseParticle *> &particles) {

}
