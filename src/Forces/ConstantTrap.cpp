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

ConstantTrap::ConstantTrap() :
				BaseForce() {
	_p_ptr = NULL;
	PBC = false;
	_r0 = -1.;
	_ref_id = -2;
}

std::tuple<std::vector<int>, std::string> ConstantTrap::init(input_file &inp) {
	BaseForce::init(inp);

	int particle;
	getInputInt(&inp, "particle", &particle, 1);
	getInputInt(&inp, "ref_particle", &_ref_id, 1);
	getInputNumber(&inp, "r0", &_r0, 1);
	getInputNumber(&inp, "stiff", &_stiff, 1);
	getInputBool(&inp, "PBC", &PBC, 0);

	int N = CONFIG_INFO->particles().size();
	if(_ref_id < 0 || _ref_id >= N) throw oxDNAException("Invalid reference particle %d for ConstantTrap", _ref_id);
	_p_ptr = CONFIG_INFO->particles()[_ref_id];

	if(particle >= N || N < -1) throw oxDNAException("Trying to add a ConstantTrap on non-existent particle %d. Aborting", particle);
	if(particle == -1) throw oxDNAException("Cannot apply ConstantTrap to all particles. Aborting");

	std::string description = Utils::sformat("ConstantTrap (stiff=%g, r0=%g, ref_particle=%d, PBC=%d", _stiff, _r0, _ref_id, PBC);

	return std::make_tuple(std::vector<int> {particle}, description);
}

LR_vector ConstantTrap::_distance(LR_vector u, LR_vector v) {
	if(PBC) {
		return CONFIG_INFO->box->min_image(u, v);
	}
	else {
		return v - u;
	}
}

LR_vector ConstantTrap::value(llint step, LR_vector &pos) {
	LR_vector dr = _distance(pos, CONFIG_INFO->box->get_abs_pos(_p_ptr)); // other - self
	number sign = copysign(1., (double) (dr.module() - _r0));
	return (_stiff * sign) * (dr / dr.module());
}

number ConstantTrap::potential(llint step, LR_vector &pos) {
	LR_vector dr = _distance(pos, CONFIG_INFO->box->get_abs_pos(_p_ptr)); // other - self
	return _stiff * fabs((dr.module() - _r0));
}
