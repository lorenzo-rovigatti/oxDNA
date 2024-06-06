/*
 * MutualTrap.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio 
 */

#include "MutualTrap.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

MutualTrap::MutualTrap() :
				BaseForce() {
	_ref_id = -2;
	_particle = -2;
	_p_ptr = NULL;
	_r0 = -1.;
	_rate = -1;
	_stiff_rate = -1;
	PBC = false;
}

std::tuple<std::vector<int>, std::string> MutualTrap::init(input_file &inp) {
	BaseForce::init(inp);

	getInputInt(&inp, "particle", &_particle, 1);
	getInputInt(&inp, "ref_particle", &_ref_id, 1);
	getInputNumber(&inp, "r0", &_r0, 1);
	getInputNumber(&inp, "stiff", &_stiff, 1);
	getInputBool(&inp, "PBC", &PBC, 0);
	_rate = 0.f; //default rate is 0
	getInputNumber(&inp, "rate", &_rate, 0);
	_stiff_rate = 0.f; //default stiff_rate is 0
	getInputNumber(&inp, "stiff_rate", &_stiff_rate, 0);

	int N = CONFIG_INFO->particles().size();
	if(_ref_id < 0 || _ref_id >= N) {
		throw oxDNAException("Invalid reference particle %d for Mutual Trap", _ref_id);
	}
	_p_ptr = CONFIG_INFO->particles()[_ref_id];

	if(_particle >= N || N < -1) {
		throw oxDNAException("Trying to add a MutualTrap on non-existent particle %d. Aborting", _particle);
	}
	if(_particle == -1) {
		throw oxDNAException("Cannot apply MutualTrap to all particles. Aborting");
	}

	std::string description = Utils::sformat("MutualTrap (stiff=%g, stiff_rate=%g, r0=%g, rate=%g, ref_particle=%d, PBC=%d)", _stiff, _stiff_rate, _r0, _rate, _ref_id, PBC);

	return std::make_tuple(std::vector<int>{_particle}, description);
}

LR_vector MutualTrap::_distance(LR_vector u, LR_vector v) {
	if(PBC) {
		return CONFIG_INFO->box->min_image(u, v);
	}
	else {
		return v - u;
	}
}

LR_vector MutualTrap::value(llint step, LR_vector &pos) {
	LR_vector dr = _distance(pos, CONFIG_INFO->box->get_abs_pos(_p_ptr));
	return (dr / dr.module()) * (dr.module() - (_r0 + (_rate * step))) * (_stiff + (_stiff_rate * step));
}

number MutualTrap::potential(llint step, LR_vector &pos) {
	LR_vector dr = _distance(pos, CONFIG_INFO->box->get_abs_pos(_p_ptr));
	return pow(dr.module() - (_r0 + (_rate * step)), 2) * ((number) 0.5) * (_stiff + (_stiff_rate * step));
}
