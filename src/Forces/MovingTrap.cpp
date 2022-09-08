/*
 * MovingTrap.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio 
 */

#include "MovingTrap.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

MovingTrap::MovingTrap() :
				BaseForce() {

}

std::tuple<std::vector<int>, std::string> MovingTrap::init(input_file &inp) {
	BaseForce::init(inp);

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	getInputNumber(&inp, "stiff", &_stiff, 1);
	getInputNumber(&inp, "rate", &_rate, 1);

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

	getInputString(&inp, "pos0", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) {
		throw oxDNAException("Could not parse pos0 %s in external forces file. Aborting", strdir.c_str());
	}
	_pos0 = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);

	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "MovingTrap");
	std::string description = Utils::sformat("MovingTrap (stiff=%g, rate=%g, dir=%g,%g,%g, pos0=%g,%g,%g", _stiff, _rate, _direction.x, _direction.y, _direction.z, _pos0.x, _pos0.y, _pos0.z);

	return std::make_tuple(particle_ids, description);
}

LR_vector MovingTrap::value(llint step, LR_vector &pos) {
	LR_vector postrap;
	number x, y, z;

	postrap.x = _pos0.x + (_rate * step) * _direction.x;
	postrap.y = _pos0.y + (_rate * step) * _direction.y;
	postrap.z = _pos0.z + (_rate * step) * _direction.z;

	x = -_stiff * (pos.x - postrap.x);
	y = -_stiff * (pos.y - postrap.y);
	z = -_stiff * (pos.z - postrap.z);

	return LR_vector(x, y, z);
}

number MovingTrap::potential(llint step, LR_vector &pos) {
	LR_vector postrap;

	postrap.x = _pos0.x + (_rate * step) * _direction.x;
	postrap.y = _pos0.y + (_rate * step) * _direction.y;
	postrap.z = _pos0.z + (_rate * step) * _direction.z;

	return (number) (0.5 * _stiff * (pos - postrap).norm());
}
