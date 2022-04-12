/*
 * LJWall.cpp
 *
 *  Created on: 04/dec/2015
 *      Author: Lorenzo
 */

#include "LJWall.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

LJWall::LJWall() :
				BaseForce() {
	_position = -1.;
	_stiff = 1;
	_n = 6;
	_sigma = 1.;
	_cutoff = 1e6;
	_generate_inside = false;
	_only_repulsive = false;
}

std::tuple<std::vector<int>, std::string> LJWall::init(input_file &inp) {
	BaseForce::init(inp);

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	getInputNumber(&inp, "stiff", &_stiff, 0);
	getInputNumber(&inp, "position", &_position, 1);
	getInputNumber(&inp, "sigma", &_sigma, 0);
	getInputInt(&inp, "n", &_n, 0);
	if(_n % 2 != 0) {
		throw oxDNAException("LJWall: n (%d) should be an even integer. Aborting", _n);
	}

	getInputBool(&inp, "only_repulsive", &_only_repulsive, 0);
	if(_only_repulsive) _cutoff = pow(2., 1. / _n);

	getInputBool(&inp, "generate_inside", &_generate_inside, 0);

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

	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "LJWall");
	std::string description = Utils::sformat("LJWall (stiff=%g, position=%g, dir=%g,%g,%g, sigma=%g, n=%d)", _stiff, _position, _direction.x, _direction.y, _direction.z, _sigma, _n);

	return std::make_tuple(particle_ids, description);
}

LR_vector LJWall::value(llint step, LR_vector &pos) {
	number distance = _direction * pos + _position; // distance from the plane
	number rel_distance = distance / _sigma; // distance from the plane in units of _sigma
	if(rel_distance > _cutoff) return LR_vector(0., 0., 0.);
	number lj_part = pow(rel_distance, -_n);
	return _direction * (4 * _n * _stiff * (2 * SQR(lj_part) - lj_part) / distance);
}

number LJWall::potential(llint step, LR_vector &pos) {
	number distance = (_direction * pos + _position) / _sigma; // distance from the plane in units of _sigma
	if(_generate_inside && distance < 0.) return 10e8;
	if(distance > _cutoff) distance = _cutoff;
	number lj_part = pow(distance, -_n);
	return 4 * _stiff * (SQR(lj_part) - lj_part);
}
