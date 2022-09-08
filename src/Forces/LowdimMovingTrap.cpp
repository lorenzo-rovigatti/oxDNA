/*
 * LowdimMovingTrap.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio 
 */

#include "LowdimMovingTrap.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

LowdimMovingTrap::LowdimMovingTrap() :
				BaseForce() {
	_visX = _visY = _visZ = false;
}

std::tuple<std::vector<int>, std::string> LowdimMovingTrap::init(input_file &inp) {
	BaseForce::init(inp);

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);
	getInputNumber(&inp, "stiff", &_stiff, 1);
	getInputNumber(&inp, "rate", &_rate, 1);

	int tmpi, tmpa[3];
	double tmpf[3];
	std::string strdir;

	getInputString(&inp, "pos0", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf, %lf, %lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) {
		throw oxDNAException("Could not parse pos0 `%s\' for LowdimMivingTrap. Aborting", strdir.c_str());
	}
	_pos0 = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);

	getInputString(&inp, "visibility", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%d, %d, %d", tmpa, tmpa + 1, tmpa + 2);
	if(tmpi != 3) {
		throw oxDNAException("Could not parse visibility `%s\' for LowdimMivingTrap. Aborting", strdir.c_str());
	}
	_visX = tmpa[0] ? true : false;
	_visY = tmpa[1] ? true : false;
	_visZ = tmpa[2] ? true : false;

	getInputString(&inp, "dir", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf, %lf, %lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) {
		throw oxDNAException("Could not parse dir `%s\' for LowdimMivingTrap. Aborting", strdir.c_str());
	}
	_direction = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	_direction.normalize();

	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "LowdimMovingTrap");
	std::string description = Utils::sformat("LowdimMovingTrap (stiff=stiffness %lf and pos=[%g,%g,%g] + (%g * t) [%g,%g,%g] and visX=%i visY=%i visZ=%i on particle %d", _stiff, _pos0.x, _pos0.y, _pos0.z, _rate, _direction.x, _direction.y, _direction.z, _visX, _visY, _visZ);

	return std::make_tuple(particle_ids, description);
}

LR_vector LowdimMovingTrap::value(llint step, LR_vector &pos) {
	LR_vector postrap;
	number x, y, z;

	postrap.x = (_visX) ? _pos0.x + (_rate * step) * _direction.x : pos.x;
	postrap.y = (_visY) ? _pos0.y + (_rate * step) * _direction.y : pos.y;
	postrap.z = (_visZ) ? _pos0.z + (_rate * step) * _direction.z : pos.z;

	x = -_stiff * (pos.x - postrap.x);
	y = -_stiff * (pos.y - postrap.y);
	z = -_stiff * (pos.z - postrap.z);

	return LR_vector(x, y, z);
}

number LowdimMovingTrap::potential(llint step, LR_vector &pos) {
	LR_vector postrap;

	postrap.x = (_visX) ? _pos0.x + (_rate * step) * _direction.x : pos.x;
	postrap.y = (_visY) ? _pos0.y + (_rate * step) * _direction.y : pos.y;
	postrap.z = (_visZ) ? _pos0.z + (_rate * step) * _direction.z : pos.z;

	return (number) (0.5 * _stiff * (pos - postrap).norm());
}
