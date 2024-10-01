/*
 * HSInteraction.cpp
 *
 *  Created on: 29/Oct/2013
 *      Author: Flavio 
 */

#include "HSInteraction.h"

HSInteraction::HSInteraction() :
				BaseInteraction() {
	ADD_INTERACTION_TO_MAP(HS, _hs_pot);
}

HSInteraction::~HSInteraction() {

}

void HSInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);
	char tmps[512];
	getInputString(&inp, "sim_type", (char*) tmps, 1);
	if(strncmp(tmps, "MC", 512) && strncmp(tmps, "MC2", 512))
		throw oxDNAException("Cannot run HS with MD");
	_rcut = (number) 1.001;
}

void HSInteraction::init() {
	_sqr_rcut = SQR(_rcut);
}

void HSInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
	for(uint i = 0; i < particles.size(); i++) {
		particles[i] = new BaseParticle();
	}
}

number HSInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, compute_r, update_forces);
}

number HSInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return (number) 0.f;
}

number HSInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	return _hs_pot(p, q, false, update_forces);
}

void HSInteraction::check_input_sanity(std::vector<BaseParticle*> &particles) {

}
