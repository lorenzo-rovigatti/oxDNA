/*
 * ConstantTrap.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio 
 */

#include "ConstantTrap.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

// constant force between two particles, to make them come together

ConstantTrap::ConstantTrap() : BaseForce() {
	_particle = -2;
	_p_ptr = NULL;
	PBC = false;
	_r0 = -1.;
	_ref_id = -2;
	_box_ptr = NULL;
}


void ConstantTrap::get_settings (input_file &inp) {
	getInputInt (&inp, "particle", &_particle, 1);
	getInputInt (&inp, "ref_particle", &_ref_id, 1);
	getInputNumber (&inp, "r0", &_r0, 1);
	getInputNumber (&inp, "stiff", &this->_stiff, 1);
	getInputBool (&inp, "PBC", &PBC, 0);
}


void ConstantTrap::init (BaseParticle ** particles, int N, BaseBox * box_ptr) {
	if (_ref_id < 0 || _ref_id >= N) throw oxDNAException ("Invalid reference particle %d for ConstantTrap", _ref_id);
	_p_ptr = particles[_ref_id];

	_box_ptr = box_ptr;
	
	if (_particle >= N || N < -1) throw oxDNAException ("Trying to add a ConstantTrap on non-existent particle %d. Aborting", _particle);
	if (_particle == -1) throw oxDNAException ("Cannot apply ConstantTrap to all particles. Aborting");

	OX_LOG (Logger::LOG_INFO, "Adding ConstantTrap (stiff=%g, r0=%g, ref_particle=%d, PBC=%d on particle %d", this->_stiff, this->_r0, _ref_id, PBC, _particle);
	particles[_particle]->add_ext_force(ForcePtr(this));
	
}

LR_vector ConstantTrap::_distance(LR_vector u, LR_vector v) {
	if (this->PBC) return _box_ptr->min_image(u, v);
	else return v - u;
}


LR_vector ConstantTrap::value (llint step, LR_vector &pos) {
	LR_vector dr = this->_distance(pos, _box_ptr->get_abs_pos(_p_ptr)); // other - self
	if (this->_site >= 0) dr -= _p_ptr->orientationT * _p_ptr->int_centers[this->_site];
	number sign = copysign (1., (double)(dr.module() - _r0));
	return (this->_stiff * sign) * (dr / dr.module());
}


number ConstantTrap::potential (llint step, LR_vector &pos) {
	LR_vector dr = this->_distance(pos, _box_ptr->get_abs_pos(_p_ptr)); // other - self
	if (this->_site >= 0) dr -= _p_ptr->orientationT * _p_ptr->int_centers[this->_site];
	return this->_stiff * fabs((dr.module () - _r0));
}

template class ConstantTrap<double>;
template class ConstantTrap<float>;
