/*
 * BoxInteraction.cpp
 *
 *  Created on: 29/Oct/2013
 *      Author: Flavio 
 */

#include "BoxInteraction.h"
#include "InteractionUtils.h"

BoxInteraction::BoxInteraction() :
				BaseInteraction() {
	ADD_INTERACTION_TO_MAP(Box, _box_pot);
}

BoxInteraction::~BoxInteraction() {

}

void BoxInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);
	char tmps[512];
	getInputString(&inp, "sim_type", (char *) tmps, 1);
	if(strncmp(tmps, "MC", 512)) throw oxDNAException("Cannot run Box with MD");

	// let's get the three dimension
	getInputString(&inp, "box_sides", (char *) tmps, 1);
	float tmpf[3];

	int tmpi = sscanf(tmps, "%g, %g, %g", tmpf, tmpf + 1, tmpf + 2);
	_sides[0] = (number) tmpf[0];
	_sides[1] = (number) tmpf[1];
	_sides[2] = (number) tmpf[2];
	if(tmpi != 3) throw oxDNAException("Could not read box sides from string \"%s\"", tmps);
	_lx = _sides[0];
	_ly = _sides[1];
	_lz = _sides[2];
	_largest = (_lx > _ly ? (_lx > _lz ? _lx : _lz) : (_ly > _lz ? _ly : _lz));
	_smallest = (_lx < _ly ? (_lx < _lz ? _lx : _lz) : (_ly < _lz ? _ly : _lz));

	OX_LOG(Logger::LOG_INFO, "Using box sides %g, %g, %g (%g, %g)", _sides[0], _sides[1], _sides[2], _smallest, _largest);

	_rcut = (number) 1.001 * sqrt(_sides[0] * _sides[0] + _sides[1] * _sides[1] + _sides[2] * _sides[2]);
	OX_LOG(Logger::LOG_INFO, "Using r_cut of %g", _rcut);
	//throw oxDNAException ("stopped here...");
}

void BoxInteraction::init() {
	_sqr_rcut = SQR(_rcut);
}

void BoxInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
	for(uint i = 0; i < particles.size(); i++) {
		particles[i] = new BaseParticle();
	}
}

number BoxInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, compute_r, update_forces);
}

number BoxInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return (number) 0.f;
}

number BoxInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	return _box_pot(p, q, false, update_forces);
}

void BoxInteraction::check_input_sanity(std::vector<BaseParticle *> &particles) {

}

bool BoxInteraction::generate_random_configuration_overlap(BaseParticle *p, BaseParticle *q) {
	LR_vector dr = _box->min_image(q->pos, p->pos);
	return InteractionUtils::box_overlap(p, q, dr, _lx, _ly, _lz);
}
