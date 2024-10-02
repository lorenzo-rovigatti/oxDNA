/*
 * RepulsiveEllipsoid.cpp
 *
 *  Created on: 26/mag/2021
 *      Author: Andrea
 */

#include "RepulsiveEllipsoid.h"

#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

RepulsiveEllipsoid::RepulsiveEllipsoid() :
				BaseForce() {
	_r_2 = LR_vector(-1., -1., -1.);
	_r_1 = LR_vector(1e-6, 1e-6, 1e-6);
	_centre = LR_vector(0., 0., 0.);
}

std::tuple<std::vector<int>, std::string> RepulsiveEllipsoid::init(input_file &inp) {
	BaseForce::init(inp);

	getInputNumber(&inp, "stiff", &_stiff, 1);

	std::string strdir;
	double tmpf[3];
	int tmpi;
	getInputString(&inp, "r_2", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) throw oxDNAException("Could not parse r0 %s in external forces file. Aborting", strdir.c_str());
	_r_2 = LR_vector  ((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);

	if(getInputString(&inp, "r_1", strdir, 0) == KEY_FOUND) {
		tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
		if(tmpi != 3) throw oxDNAException("Could not parse r_ext %s in external forces file. Aborting", strdir.c_str());
		_r_1 = LR_vector  ((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	}

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	if(getInputString(&inp, "center", strdir, 0) == KEY_FOUND) {
		tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
		if(tmpi != 3) throw oxDNAException("Could not parse center %s in external forces file. Aborting", strdir.c_str());
		_centre = LR_vector  ((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	}

	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "RepulsiveEllipse");
	std::string description = Utils::sformat("RepulsiveEllipse(stiff=%g, r0=%g,%g,%g, center=%g,%g,%g)", _stiff, _r_2.x, _r_2.y, _r_2.z, _centre.x, _centre.y, _centre.z);

	return std::make_tuple(particle_ids, description);
}

LR_vector RepulsiveEllipsoid::value(llint step, LR_vector &pos) {
	LR_vector dist = CONFIG_INFO->box->min_image(_centre, pos);

	double internal_cut = SQR(dist.x) / SQR(_r_2.x) + SQR(dist.y) / SQR(_r_2.y) + SQR(dist.z) / SQR(_r_2.z);
	double external_cut = SQR(dist.x) / SQR(_r_1.x) + SQR(dist.y) / SQR(_r_1.y) + SQR(dist.z) / SQR(_r_1.z);

	if(internal_cut < 1. && external_cut > 1.) return LR_vector(0., 0., 0.);
	else {
		dist.normalize();
		return -_stiff * dist;
	}
}

number RepulsiveEllipsoid::potential(llint step, LR_vector &pos) {
	LR_vector dist = CONFIG_INFO->box->min_image(_centre, pos);

	double internal_cut = SQR(dist.x) / SQR(_r_2.x) + SQR(dist.y) / SQR(_r_2.y) + SQR(dist.z) / SQR(_r_2.z);
	double external_cut = SQR(dist.x) / SQR(_r_1.x) + SQR(dist.y) / SQR(_r_1.y) + SQR(dist.z) / SQR(_r_1.z);
	if(internal_cut < 1. && external_cut > 1.) return 0.;
	else {
		return _stiff * dist.module();
	}
}
