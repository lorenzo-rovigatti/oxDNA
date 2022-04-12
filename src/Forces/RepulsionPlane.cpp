/*
 * RepulsionPlane.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio 
 */

#include "RepulsionPlane.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

RepulsionPlane::RepulsionPlane() :
				BaseForce() {
	_position = -1.;
}

std::tuple<std::vector<int>, std::string> RepulsionPlane::init(input_file &inp) {
	BaseForce::init(inp);

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	getInputNumber(&inp, "stiff", &_stiff, 1);
	getInputNumber(&inp, "position", &_position, 1);

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

	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "RepulsionPlane");
	std::string description = Utils::sformat("RepulsionPlane (stiff=%g, position=%g, dir=%g,%g,%g", _stiff, _position, _direction.x, _direction.y, _direction.z);

	return std::make_tuple(particle_ids, description);
}

LR_vector RepulsionPlane::value(llint step, LR_vector &pos) {
	number distance = _direction * pos + _position;
	if(distance >= 0.) return LR_vector(0., 0., 0.);
	else return -(distance * _stiff) * _direction;
}

number RepulsionPlane::potential(llint step, LR_vector &pos) {
	number distance = _direction * pos + _position; // distance from the plane
	if(distance >= 0.) return 0.;
	else return (number) (0.5 * _stiff * SQR(distance));
}
