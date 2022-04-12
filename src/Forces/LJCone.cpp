/*
 * LJCone.cpp
 *
 *  Created on: 24/aug/2017
 *      Author: Lorenzo
 */

#include "LJCone.h"

#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

LJCone::LJCone() :
				BaseForce() {
	_stiff = 1;
	_n = 6;
	_sigma = 1.;
	_cutoff = 1e6;
	_alpha = _sin_alpha = _cos_alpha = _tan_alpha = 0.;
	_generate_inside = false;
	_only_repulsive = false;
}

std::tuple<std::vector<int>, std::string> LJCone::init(input_file &inp) {
	BaseForce::init(inp);

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	getInputNumber(&inp, "stiff", &_stiff, 0);
	getInputNumber(&inp, "sigma", &_sigma, 0);
	getInputNumber(&inp, "alpha", &_alpha, 1);
	getInputInt(&inp, "n", &_n, 0);
	if(_n % 2 != 0) throw oxDNAException("RepulsiveCone: n (%d) should be an even integer. Aborting", _n);

	getInputBool(&inp, "only_repulsive", &_only_repulsive, 0);
	if(_only_repulsive) _cutoff = pow(2., 1. / _n);

	getInputBool(&inp, "generate_inside", &_generate_inside, 0);

	int tmpi;
	double tmpf[3];
	std::string str_vector;
	getInputString(&inp, "dir", str_vector, 1);
	tmpi = sscanf(str_vector.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) {
		throw oxDNAException("RepulsiveCone: could not parse dir %s in external forces file. Aborting", str_vector.c_str());
	}
	_direction = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	_direction.normalize();

	getInputString(&inp, "pos0", str_vector, 1);
	tmpi = sscanf(str_vector.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) {
		throw oxDNAException("RepulsiveCone: could not parse pos0 %s in external forces file. Aborting", str_vector.c_str());
	}
	_pos0 = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);

	_sin_alpha = sin(_alpha);
	_cos_alpha = cos(_alpha);
	_tan_alpha = tan(_alpha);

	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "ConstantRateForce");
	std::string description = Utils::sformat("RepulsiveCone (stiff=%g, pos0=%g,%g,%g, dir=%g,%g,%g, sigma=%g, n=%d)", _stiff, _pos0.x, _pos0.y, _pos0.z, _direction.x, _direction.y, _direction.z, _sigma, _n);

	return std::make_tuple(particle_ids, description);
}

LR_vector LJCone::value(llint step, LR_vector &pos) {
	LR_vector v_from_apex = pos - _pos0;

	number d_along_axis = v_from_apex * _direction;
	LR_vector v_along_axis = _direction * d_along_axis;

	LR_vector v_from_axis = v_along_axis - v_from_apex;
	number d_from_axis = v_from_axis.module();
	number d_from_cone = d_along_axis * _sin_alpha - d_from_axis * _cos_alpha;
	number rel_distance = d_from_cone / _sigma;

	if(rel_distance > _cutoff) return LR_vector(0., 0., 0.);

	// now we compute the normal to the cone
	number C = d_from_axis * _tan_alpha;
	LR_vector C_v = (d_along_axis + C) * _direction;
	LR_vector normal = C_v - v_from_apex;
	normal.normalize();

	number lj_part = pow(rel_distance, -_n);
	return normal * (4 * _n * _stiff * (2 * SQR(lj_part) - lj_part) / d_from_cone);
}

number LJCone::potential(llint step, LR_vector &pos) {
	LR_vector v_from_apex = pos - _pos0;

	number d_along_axis = v_from_apex * _direction;
	LR_vector v_along_axis = _direction * d_along_axis;
	LR_vector v_from_axis = v_from_apex - v_along_axis;
	number d_from_cone = d_along_axis * _sin_alpha - v_from_axis.module() * _cos_alpha;
	number rel_distance = d_from_cone / _sigma; // distance from the plane in units of _sigma

	if(_generate_inside && rel_distance < 0.) return 10e8;
	if(rel_distance > _cutoff) rel_distance = _cutoff;
	number lj_part = pow(rel_distance, -_n);
	return 4 * _stiff * (SQR(lj_part) - lj_part);
}
