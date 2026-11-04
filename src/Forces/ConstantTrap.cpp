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

	std::string particle_string;
	std::string ref_particle_string;
	getInputString(&inp, "particle", particle_string, 1);
	getInputString(&inp, "ref_particle", ref_particle_string, 1);
	int particle = Utils::get_single_particle_from_string(CONFIG_INFO->particles(), particle_string, "ConstantTrap particle");
	_ref_id = Utils::get_single_particle_from_string(CONFIG_INFO->particles(), ref_particle_string, "ConstantTrap ref_particle");
	getInputNumber(&inp, "r0", &_r0, 1);
	getInputNumber(&inp, "stiff", &_stiff, 1);
	getInputBool(&inp, "PBC", &PBC, 0);

	_p_ptr = CONFIG_INFO->particles()[_ref_id];

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

LR_vector ConstantTrap::force(llint step, LR_vector &pos) {
	LR_vector dr = _distance(pos, CONFIG_INFO->box->get_abs_pos(_p_ptr)); // other - self
	number sign = copysign(1., (double) (dr.module() - _r0));
	return (_stiff * sign) * (dr / dr.module());
}

number ConstantTrap::potential(llint step, LR_vector &pos) {
	LR_vector dr = _distance(pos, CONFIG_INFO->box->get_abs_pos(_p_ptr)); // other - self
	return _stiff * fabs((dr.module() - _r0));
}
