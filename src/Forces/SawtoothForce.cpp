/*
 * ExternalForce.cpp
 *
 *  Created on: 8/apr/2015
 *      Author: Megan
 */

#include "SawtoothForce.h"
#include "../Particles/BaseParticle.h"

SawtoothForce::SawtoothForce() :
				BaseForce() {
	_wait_time = 0.0;
	_increment = 0.0;
}

SawtoothForce::~SawtoothForce() {

}

std::tuple<std::vector<int>, std::string> SawtoothForce::init(input_file &inp) {
	BaseForce::init(inp);

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	getInputNumber(&inp, "F0", &_F0, 1);
	getInputNumber(&inp, "wait_time", &_wait_time, 1);
	getInputNumber(&inp, "increment", &_increment, 1);

	std::string strdir;
	getInputString(&inp, "dir", strdir, 1);
	double x, y, z;
	int tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", &x, &y, &z);
	if(tmpi != 3) throw oxDNAException("could not parse direction in external_forces_file. Dying badly");
	_direction = LR_vector((number) x, (number) y, number(z));
	_direction.normalize();

	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "SawtoothForce");
	std::string description = Utils::sformat("SawtoothForce (F==%g, wait_time=%g, increment=%g, dir=%g,%g,%g", _F0, _wait_time, _increment, _direction.x, _direction.y, _direction.z);

	return std::make_tuple(particle_ids, description);
}

LR_vector SawtoothForce::value(llint step, LR_vector &pos) {
	number x = (_F0 + (int) ((step - 1) / _wait_time) * _increment) * _direction.x;
	number y = (_F0 + (int) ((step - 1) / _wait_time) * _increment) * _direction.y;
	number z = (_F0 + (int) ((step - 1) / _wait_time) * _increment) * _direction.z;
	return LR_vector(x, y, z);
}

number SawtoothForce::potential(llint step, LR_vector &pos) {
	return (number) -(_F0 + (int) ((step - 1) / _wait_time) * _increment) * (pos * _direction);
}
