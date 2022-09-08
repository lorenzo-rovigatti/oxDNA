/*
 * HardWall.cpp
 *
 *  Created on: 5/oct/2016
 *      Author: Lorenzo
 */

#include "HardWall.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

HardWall::HardWall() :
				BaseForce() {
	_position = -1.;
	_stiff = 1;
	_sigma = 1.;
}

std::tuple<std::vector<int>, std::string> HardWall::init(input_file &inp) {
	BaseForce::init(inp);

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	getInputNumber(&inp, "position", &_position, 1);
	getInputNumber(&inp, "sigma", &_sigma, 0);

	int tmpi;
	double tmpf[3];
	std::string strdir;
	getInputString(&inp, "dir", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) {
		throw oxDNAException("Could not parse dir %s in external forces file. Aborting", strdir.c_str());
	}
	_direction = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	_direction.normalize();

	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "HardWall");
	std::string description = Utils::sformat("HardWall (position=%g, dir=%g,%g,%g, sigma=%g)", _position, _direction.x, _direction.y, _direction.z, _sigma);

	return std::make_tuple(particle_ids, description);
}

LR_vector HardWall::value(llint step, LR_vector &pos) {
	throw oxDNAException("HardWall can be used only in Monte Carlo simulations.");
}

number HardWall::potential(llint step, LR_vector &pos) {
	number distance = (_direction * pos + _position) / _sigma; // distance from the plane in units of _sigma
	if(distance < 1.) return 10e8;
	return 0.;
}
