/*
 * RepulsiveSphere.cpp
 *
 *  Created on: 28/nov/2014
 *      Author: Flavio 
 */

#include "RepulsiveSphere.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

RepulsiveSphere::RepulsiveSphere() :
				BaseForce() {
	_r0 = -1.;
	_r_ext = 1e10;
	_center = LR_vector(0., 0., 0.);
	_rate = 0.;
}

std::tuple<std::vector<int>, std::string> RepulsiveSphere::init(input_file &inp) {
	BaseForce::init(inp);

	getInputNumber(&inp, "stiff", &_stiff, 1);
	getInputNumber(&inp, "r0", &_r0, 1);
	getInputNumber(&inp, "rate", &_rate, 0);
	getInputNumber(&inp, "r_ext", &_r_ext, 0);

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	std::string strdir;
	if(getInputString(&inp, "center", strdir, 0) == KEY_FOUND) {
		double tmpf[3];
		int tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
		if(tmpi != 3) throw oxDNAException("Could not parse center %s in external forces file. Aborting", strdir.c_str());
		_center = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	}

	std::string description = Utils::sformat("RepulsiveSphere (stiff=%g, r0=%g, rate=%g, center=%g,%g,%g)", _stiff, _r0, _rate, _center.x, _center.y, _center.z);
	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "RepulsiveSphere");

	return std::make_tuple(particle_ids, description);
}

LR_vector RepulsiveSphere::value(llint step, LR_vector &pos) {
	LR_vector dist = CONFIG_INFO->box->min_image(_center, pos);
	number mdist = dist.module();
	number radius = _r0 + _rate * (number) step;

	if(mdist <= radius || mdist >= _r_ext) return LR_vector(0., 0., 0.);
	else return dist * (-_stiff * (1. - radius / mdist));
}

number RepulsiveSphere::potential(llint step, LR_vector &pos) {
	LR_vector dist = CONFIG_INFO->box->min_image(_center, pos);
	number mdist = dist.module();
	number radius = _r0 + _rate * (number) step;

	if(mdist <= radius || mdist >= _r_ext) return 0.;
	else return 0.5 * _stiff * (mdist - radius) * (mdist - radius);
}
