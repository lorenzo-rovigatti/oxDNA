/*
 * BoxInteraction.cpp
 *
 *  Created on: 29/Oct/2013
 *      Author: Flavio 
 */

#include "BoxInteraction.h"
#include "InteractionUtils.h"


BoxInteraction::BoxInteraction() : BaseInteraction<number, BoxInteraction >() {
	this->_int_map[Box] = &BoxInteraction::_box_pot;
}


BoxInteraction::~BoxInteraction() {

}


void BoxInteraction::get_settings(input_file &inp) {
	IBaseInteraction::get_settings(inp);
	char tmps[512];
	getInputString (&inp, "sim_type", (char *)tmps, 1);
	if (strncmp(tmps, "MC", 512)) throw oxDNAException ("Cannot run Box with MD");
	
	// let's get the three dimension
	getInputString (&inp, "box_sides", (char *)tmps, 1);
	float tmpf[3];
	
	int tmpi = sscanf (tmps, "%g, %g, %g", tmpf, tmpf + 1, tmpf + 2);
	_sides[0] = (number) tmpf[0];
	_sides[1] = (number) tmpf[1];
	_sides[2] = (number) tmpf[2];
	if (tmpi != 3) throw oxDNAException ("Could not read box sides from string \"%s\"", tmps);
	_lx = _sides[0];
	_ly = _sides[1];
	_lz = _sides[2];
	_largest = (_lx > _ly ? (_lx > _lz ? _lx : _lz) : (_ly > _lz ? _ly : _lz));
	_smallest = (_lx < _ly ? (_lx < _lz ? _lx : _lz) : (_ly < _lz ? _ly : _lz));
	
	OX_LOG(Logger::LOG_INFO, "Using box sides %g, %g, %g (%g, %g)", _sides[0], _sides[1], _sides[2], _smallest, _largest);
	
	this->_rcut = (number) 1.001 * sqrt(_sides[0] * _sides[0] + _sides[1] * _sides[1] + _sides[2] * _sides[2]);
	OX_LOG(Logger::LOG_INFO, "Using r_cut of %g", this->_rcut);
	//throw oxDNAException ("stopped here...");
}


void BoxInteraction::init() {
	this->_sqr_rcut = SQR(this->_rcut);
}


void BoxInteraction::allocate_particles(BaseParticle **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new BaseParticle();
}


number BoxInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}


number BoxInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {
	return (number) 0.f;
}


number BoxInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces) {
	LR_vector computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}
	
	return _box_pot (p, q, r, update_forces);
}


void BoxInteraction::check_input_sanity(BaseParticle **particles, int N) {

}


bool BoxInteraction::generate_random_configuration_overlap (BaseParticle *p, BaseParticle *q) {
	LR_vector dr = this->_box->min_image (q->pos, p->pos);
	return InteractionUtils::box_overlap (p, q, dr, _lx, _ly, _lz);
}

template class BoxInteraction<float>;
template class BoxInteraction<double>;
