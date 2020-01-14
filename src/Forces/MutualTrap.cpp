/*
 * MutualTrap.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio 
 */

#include "MutualTrap.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"


MutualTrap::MutualTrap() : BaseForce() {
	_ref_id = -2;
	_particle = -2;
	_p_ptr = NULL;
	_r0 = -1.;
	PBC = false;
	_box_ptr = NULL;
}


void MutualTrap::get_settings (input_file &inp) {
	getInputInt (&inp, "particle", &_particle, 1);
	getInputInt (&inp, "ref_particle", &_ref_id, 1);
	getInputNumber (&inp, "r0", &_r0, 1);
	getInputNumber (&inp, "stiff", &this->_stiff, 1);
	getInputBool (&inp, "PBC", &PBC, 0);
	this->_rate = 0.f; //default rate is 0
	getInputNumber(&inp, "rate", &this->_rate, 0);
}


void MutualTrap::init (std::vector<BaseParticle *> & particles, int N, BaseBox * box_ptr){
	if (_ref_id < 0 || _ref_id >= N) throw oxDNAException ("Invalid reference particle %d for Mutual Trap", _ref_id);
	_p_ptr = particles[_ref_id];

	this->_box_ptr = box_ptr;
	
	if(_particle >= N || N < -1) throw oxDNAException ("Trying to add a MutualTrap on non-existent particle %d. Aborting", _particle);
	if(_particle == -1) throw oxDNAException ("Cannot apply MutualTrap to all particles. Aborting");

	OX_LOG (Logger::LOG_INFO, "Adding MutualTrap (stiff=%g, rate=%g, r0=%g, ref_particle=%d, PBC=%d on particle %d", this->_stiff, this->_rate, this->_r0, _ref_id, PBC, _particle);
	particles[_particle]->add_ext_force(ForcePtr(this));
	
}


LR_vector MutualTrap::_distance(LR_vector u, LR_vector v) {
	if (this->PBC) return this->_box_ptr->min_image(u, v);
	else return v - u;
}


LR_vector MutualTrap::value (llint step, LR_vector &pos) {
	LR_vector dr = this->_distance(pos, this->_box_ptr->get_abs_pos(_p_ptr));
	return (dr / dr.module()) * (dr.module() - (_r0 + (this->_rate * step))) * this->_stiff;
}


number MutualTrap::potential (llint step, LR_vector &pos) {
	LR_vector dr = this->_distance(pos, this->_box_ptr->get_abs_pos(_p_ptr));
	return pow (dr.module() - (_r0 + (this->_rate * step)), 2) * ((number) 0.5) * this->_stiff;
}
