/*
 * RepulsiveEllipse.cpp
 *
 *  Created on: 26/mag/2021
 *      Author: Andrea
 */

#include "RepulsiveEllipse.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

RepulsiveEllipse::RepulsiveEllipse() :
				BaseForce() {
	_r_2 = LR_vector(-1., -1., -1.);
	_r_1 = LR_vector(1e-6, 1e-6, 1e-6);
	_center = LR_vector(0., 0., 0.);
	_box_ptr = nullptr;
}

std::tuple<std::vector<int>, std::string> RepulsiveEllipse::init(input_file &inp, BaseBox *box_ptr) {
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
		_center = LR_vector  ((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	}

	_box_ptr = box_ptr;

	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "RepulsiveEllipse");
	std::string description = Utils::sformat("RepulsiveEllipse(stiff=%g, r0=%g,%g,%g, center=%g,%g,%g)", _stiff, _r_2.x, _r_2.y, _r_2.z, _center.x, _center.y, _center.z);

	return std::make_tuple(particle_ids, description);
}

LR_vector RepulsiveEllipse::value(llint step, LR_vector &pos) {
	LR_vector dist = _box_ptr->min_image(_center, pos);
	LR_vector radius = _r_2;

	double internal_cut = SQR(dist.x) / SQR(radius.x) + SQR(dist.y) / SQR(radius.y) + SQR(dist.z) / SQR(radius.z);
	double external_cut = SQR(dist.x) / SQR(_r_1.x) + SQR(dist.y) / SQR(_r_1.y) + SQR(dist.z) / SQR(_r_1.z);

	if(internal_cut < 1. && external_cut > 1.) return LR_vector(0., 0., 0.);
	else {
		dist.normalize();
		return -_stiff * dist;
	}
}

number RepulsiveEllipse::potential(llint step, LR_vector &pos) {
	LR_vector dist = _box_ptr->min_image(_center, pos);
	LR_vector radius = _r_2;

	double internal_cut = SQR(dist.x) / SQR(radius.x) + SQR(dist.y) / SQR(radius.y) + SQR(dist.z) / SQR(radius.z);
	double external_cut = SQR(dist.x) / SQR(_r_1.x) + SQR(dist.y) / SQR(_r_1.y) + SQR(dist.z) / SQR(_r_1.z);
	if(internal_cut < 1. && external_cut > 1.) return 0.;
	else {
		return _stiff * dist.module();
	}
}
