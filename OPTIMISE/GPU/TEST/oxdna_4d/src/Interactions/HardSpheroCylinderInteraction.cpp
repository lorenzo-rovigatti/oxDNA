/*
 * HardSpheroCylinderInteraction.cpp
 *
 *  Created on: 29/Oct/2013
 *      Author: Flavio 
 */

#include "HardSpheroCylinderInteraction.h"

HardSpheroCylinderInteraction::HardSpheroCylinderInteraction() :
				BaseInteraction() {
	ADD_INTERACTION_TO_MAP(0, _hsc_pot);
}

HardSpheroCylinderInteraction::~HardSpheroCylinderInteraction() {

}

void HardSpheroCylinderInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);
	char tmps[512];
	getInputString(&inp, "sim_type", (char*) tmps, 1);
	if(strncmp(tmps, "MC", 512))
		throw oxDNAException("Cannot run HardSpheroCylinder with MD");

	// length
	float tmpf;
	getInputFloat(&inp, "length", &tmpf, 1);
	_length = (number) tmpf;
	if(_length < 0.)
		throw oxDNAException("Cannot run hard spherocylinders with negative lenght");

	// rcut is 2. * (radius + length)
	_rcut = (number) 1.001 + _length;

	OX_LOG(Logger::LOG_INFO, "Initializing HardSpheroCylinder interaction with length %g [_rcut = %g]", _length, _rcut);
}

void HardSpheroCylinderInteraction::init() {
	_sqr_rcut = SQR(_rcut);
}

void HardSpheroCylinderInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
	for(uint i = 0; i < particles.size(); i++) {
		particles[i] = new BaseParticle();
	}
}

number HardSpheroCylinderInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, compute_r, update_forces);
}

number HardSpheroCylinderInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	return (number) 0.f;
}

number HardSpheroCylinderInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	return _hsc_pot(p, q, false, update_forces);
}

void HardSpheroCylinderInteraction::check_input_sanity(std::vector<BaseParticle*> &particles) {

}
