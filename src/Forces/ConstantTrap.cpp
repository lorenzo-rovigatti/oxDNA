/*
 * ConstantTrap.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio 
 */

#include "ConstantTrap.h"
#include "../Particles/BaseParticle.h"

// constant force between two particles, to make them come together
template<typename number>
ConstantTrap<number>::ConstantTrap() : BaseForce<number>() {
	_particle = -2;
	_p_ptr = NULL;
	PBC = false;
	_r0 = -1.;
	_ref_id = -2;
}

template <typename number>
void ConstantTrap<number>::get_settings (input_file &inp) {
	getInputInt (&inp, "particle", &_particle, 1);
	getInputInt (&inp, "ref_particle", &_ref_id, 1);
	getInputNumber (&inp, "r0", &_r0, 1);
	getInputNumber (&inp, "stiff", &this->_stiff, 1);
	getInputBool (&inp, "PBC", &PBC, 0);
}

template <typename number>
void ConstantTrap<number>::init (BaseParticle<number> ** particles, int N, number * my_box_side_ptr){
	if (_ref_id < 0 || _ref_id >= N) throw oxDNAException ("Invalid reference particle %d for ConstantTrap", _ref_id);
	_p_ptr = particles[_ref_id];

	this->box_side_ptr = my_box_side_ptr;
	
	if (_particle >= N || N < -1) throw oxDNAException ("Trying to add a ConstantTrap on non-existent particle %d. Aborting", _particle);
	if (_particle == -1) throw oxDNAException ("Cannot apply ConstantTrap to all particles. Aborting");

	OX_LOG (Logger::LOG_INFO, "Adding ConstantTrap (stiff=%g, r0=%g, ref_particle=%d, PBC=%d on particle %d", this->_stiff, this->_r0, _ref_id, PBC, _particle);
	particles[_particle]->add_ext_force(this);
	
}
template<typename number>
LR_vector<number> ConstantTrap<number>::_distance(LR_vector<number> u, LR_vector<number> v) {
	if (this->PBC) return v.minimum_image(u, *(this->box_side_ptr));
	else return v - u;
}

template<typename number>
LR_vector<number> ConstantTrap<number>::value (llint step, LR_vector<number> &pos) {
	LR_vector<number> dr = this->_distance(pos, _p_ptr->get_abs_pos(*(this->box_side_ptr))); // other - self
	if (this->_site >= 0) dr -= _p_ptr->orientationT * _p_ptr->int_centers[this->_site];
	number sign = copysign (1., (double)(dr.module() - _r0));
	return (this->_stiff * sign) * (dr / dr.module());
}

template <typename number>
number ConstantTrap<number>::potential (llint step, LR_vector<number> &pos) {
	LR_vector<number> dr = this->_distance(pos, _p_ptr->get_abs_pos(*(this->box_side_ptr))); // other - self
	if (this->_site >= 0) dr -= _p_ptr->orientationT * _p_ptr->int_centers[this->_site];
	return this->_stiff * fabs((dr.module () - _r0));
}

template class ConstantTrap<double>;
template class ConstantTrap<float>;
