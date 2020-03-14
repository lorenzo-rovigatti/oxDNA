/*
 * HardSpheroCylinderInteraction.cpp
 *
 *  Created on: 29/Oct/2013
 *      Author: Flavio 
 */

#include "HardSpheroCylinderInteraction.h"

HardSpheroCylinderInteraction::HardSpheroCylinderInteraction() :
				BaseInteraction<HardSpheroCylinderInteraction>() {
	this->_int_map[0] = &HardSpheroCylinderInteraction::_hsc_pot;
}

HardSpheroCylinderInteraction::~HardSpheroCylinderInteraction() {

}

void HardSpheroCylinderInteraction::get_settings(input_file &inp) {
	IBaseInteraction::get_settings(inp);
	char tmps[512];
	getInputString(&inp, "sim_type", (char *) tmps, 1);
	if(strncmp(tmps, "MC", 512)) throw oxDNAException("Cannot run HardSpheroCylinder with MD");

	// length
	float tmpf;
	getInputFloat(&inp, "length", &tmpf, 1);
	_length = (number) tmpf;
	if(_length < 0.) throw oxDNAException("Cannot run hard spherocylinders with negative lenght");

	// rcut is 2. * (radius + length)
	this->_rcut = (number) 1.001 + _length;

	OX_LOG(Logger::LOG_INFO, "Initializing HardSpheroCylinder interaction with length %g [_rcut = %g]", _length, this->_rcut);
}

void HardSpheroCylinderInteraction::init() {
	this->_sqr_rcut = SQR(this->_rcut);
}

void HardSpheroCylinderInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
	for(uint i = 0; i < particles.size(); i++) {
		particles[i] = new BaseParticle();
	}
}

number HardSpheroCylinderInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}

number HardSpheroCylinderInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {
	return (number) 0.f;
}

number HardSpheroCylinderInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {
	LR_vector computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	return _hsc_pot(p, q, r, update_forces);
}

void HardSpheroCylinderInteraction::check_input_sanity(std::vector<BaseParticle *> &particles) {

}
