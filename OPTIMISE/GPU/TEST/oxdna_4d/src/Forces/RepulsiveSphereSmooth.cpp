/*
 * RepulsiveSphere.cpp
 *
 *  Created on: 28/nov/2014
 *      Author: Flavio 
 */

#include "RepulsiveSphereSmooth.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

RepulsiveSphereSmooth::RepulsiveSphereSmooth() :
				BaseForce() {
	_r0 = -1.;
	_r_ext = 1e10;
	_center = LR_vector(0., 0., 0.);
	_alpha = -1.;
	_smooth = -1.;
}

std::tuple<std::vector<int>, std::string> RepulsiveSphereSmooth::init(input_file &inp) {
	BaseForce::init(inp);

	getInputNumber(&inp, "stiff", &_stiff, 1);
	getInputNumber(&inp, "r0", &_r0, 1);
	getInputNumber(&inp, "r_ext", &_r_ext, 0);
	getInputNumber(&inp, "r_ext", &_smooth, 1);
	getInputNumber(&inp, "r_ext", &_alpha, 1);

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	std::string strdir;
	if(getInputString(&inp, "center", strdir, 0) == KEY_FOUND) {
		double tmpf[3];
		int tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
		if(tmpi != 3) throw oxDNAException("Could not parse center %s in external forces file. Aborting", strdir.c_str());
		_center = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	}

	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "RepulsiveSphereSmooth");
	std::string description = Utils::sformat("RepulsiveSphereSmooth force (stiff=%g, r0=%g,  center=%g,%g,%g)", _stiff, _r0, _center.x, _center.y, _center.z);

	return std::make_tuple(particle_ids, description);
}

LR_vector RepulsiveSphereSmooth::value(llint step, LR_vector &pos) {
	LR_vector dist = CONFIG_INFO->box->min_image(_center, pos);
	number mdist = dist.module();

	if(mdist < _r0 || mdist > _r_ext) {
		return LR_vector(0., 0., 0.);
	}
	else {
		if(mdist >= _alpha && mdist <= _r_ext) {
			return dist * (-(_stiff * 0.5 * exp((mdist - _alpha) / _smooth)) / mdist);
		}
		else {
			return dist * (-(_stiff * mdist - _stiff * 0.5 * exp(-(mdist - _alpha) / _smooth)) / mdist);
		}
	}
}

number RepulsiveSphereSmooth::potential(llint step, LR_vector &pos) {
	LR_vector dist = CONFIG_INFO->box->min_image(_center, pos);
	number mdist = dist.module();

	if(mdist < _r0 || mdist > _r_ext) return 0.;
	else {
		if(mdist >= _alpha && mdist <= _r_ext) return (_stiff * 0.5 * _smooth * exp((mdist - _alpha) / _smooth));
		else return (_stiff * mdist + _stiff * 0.5 * _smooth * exp(-(mdist - _alpha) / _smooth));
	}
}
