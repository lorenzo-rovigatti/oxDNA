/*
 * ConstantRateForce.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio 
 *              Modified by Ferdinando on 11/Dec/2015
 */

#include "ConstantRateForce.h"
#include "../Particles/BaseParticle.h"

using namespace std;

ConstantRateForce::ConstantRateForce() :
				BaseForce() {
	dir_as_centre = false;
}

ConstantRateForce::~ConstantRateForce() {

}

std::tuple<std::vector<int>, std::string> ConstantRateForce::init(input_file &inp) {
	BaseForce::init(inp);

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	getInputNumber(&inp, "F0", &_F0, 1);
	getInputNumber(&inp, "rate", &_rate, 1);
	getInputBool(&inp, "dir_as_centre", &dir_as_centre, 0);

	string strdir;
	getInputString(&inp, "dir", strdir, 1);
	vector<string> spl = Utils::split(strdir, ',');
	if(spl.size() != 3) {
		throw oxDNAException("Could not parse 'dir' in external_forces_file. Dying badly");
	}

	_direction.x = atof(spl[0].c_str());
	_direction.y = atof(spl[1].c_str());
	_direction.z = atof(spl[2].c_str());

	if(!dir_as_centre) {
		_direction.normalize();
	}

	std::string description = Utils::sformat("ConstantRateForce (F=%g, rate=%g, dir=%g,%g,%g)", _F0, _rate, _direction.x, _direction.y, _direction.z);
	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "ConstantRateForce");

	return std::make_tuple(particle_ids, description);
}

LR_vector ConstantRateForce::value(llint step, LR_vector &pos) {
	LR_vector dir = _direction;
	if(dir_as_centre) {
		dir -= pos;
		dir.normalize();
	}
	LR_vector force = (_F0 + _rate * step) * dir;
	return force;
}

number ConstantRateForce::potential(llint step, LR_vector &pos) {
	number strength = -(_F0 + _rate * step);
	if(dir_as_centre) {
		number dist = (_direction - pos).module();
		return strength * dist;
	}
	return strength * (pos * _direction);
}
