/*
 * YukawaSphere.cpp
 *
 *  Created on: 31/may/2024
 *      Author: Lorenzo
 */

#include "YukawaSphere.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

YukawaSphere::YukawaSphere() :
				BaseForce() {

}

std::tuple<std::vector<int>, std::string> YukawaSphere::init(input_file &inp) {
	BaseForce::init(inp);

	getInputNumber(&inp, "radius", &_radius, 1);
	getInputNumber(&inp, "WCA_epsilon", &_epsilon, 0);
	getInputNumber(&inp, "WCA_sigma", &_sigma, 0);
	getInputNumber(&inp, "WCA_n", &_sigma, 0);
	_WCA_cutoff = _sigma * pow(2., 1. / _WCA_n);

	getInputNumber(&inp, "debye_length", &_debye_length, 1);
	getInputNumber(&inp, "debye_A", &_debye_A, 1);
	if(getInputNumber(&inp, "cutoff", &_cutoff, 0) == KEY_NOT_FOUND) {
		_cutoff = 4.0 * _debye_length;
	}

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	std::string strdir;
	if(getInputString(&inp, "center", strdir, 0) == KEY_FOUND) {
		auto center_as_vector = Utils::split_to_numbers(strdir, ",");
		if(center_as_vector.size() == 3) {
			_center = LR_vector(center_as_vector[0], center_as_vector[1], center_as_vector[2]);
		}
		else {
			throw oxDNAException("Could not parse center %s in external forces file. Aborting", strdir.c_str());
		}
	}

	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "YukawaSphere");
	std::string description = Utils::sformat("YukawaSphere force (radius=%g, center=%g,%g,%g, eps=%g, sigma=%g, debye_length=%g, debye_A=%g, cutoff=%g)", _radius, _center.x, _center.y, _center.z, _epsilon, _sigma, _debye_length, _debye_A, _cutoff);

	return std::make_tuple(particle_ids, description);
}

LR_vector YukawaSphere::value(llint step, LR_vector &pos) {
	LR_vector dist = CONFIG_INFO->box->min_image(_center, pos);
	number dist_surface = _radius - dist.module();
	if(dist_surface < 0) {
		throw oxDNAException("Found a particle beyond the Yukawa sphere radius %lf. Check your initial configuration or your timestep.", _radius);
	}

	LR_vector force;
	if(dist_surface < _cutoff) {
		number dist_surface_sqr = SQR(dist_surface);
		LR_vector direction = -dist / dist.module();

		force = direction * ((_debye_A * exp(-dist_surface / _debye_length)) * (1.0 / (dist_surface * _debye_length) + 1.0 / (dist_surface_sqr)));

		if(dist_surface < _WCA_cutoff) {
			number WCA_part = pow(_sigma / dist_surface, 6);
			force += direction * (4 * _epsilon * _WCA_n * (2 * SQR(WCA_part) - WCA_part) / dist_surface);
		}
	}

	return force;
}

number YukawaSphere::potential(llint step, LR_vector &pos) {
	number sqr_r = CONFIG_INFO->box->sqr_min_image_distance(_center, pos);
	number dist_surface = _radius - sqrt(sqr_r);
	if(dist_surface < 0) {
		return 1e10;
	}

	number energy = 0.;
	if(dist_surface < _cutoff) {
		energy = exp(-dist_surface / _debye_length) * _debye_A / dist_surface;
		if(dist_surface < _WCA_cutoff) {
			number WCA_part = pow(_sigma / dist_surface, 6);
			energy += 4 * _epsilon * (SQR(WCA_part) - WCA_part) + _epsilon;

		}
	}
	return energy;
}
