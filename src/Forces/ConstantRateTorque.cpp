/*
 * ConstantRateTorque.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio 
 */

#include "ConstantRateTorque.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

ConstantRateTorque::ConstantRateTorque() :
				BaseForce() {

}

std::tuple<std::vector<int>, std::string> ConstantRateTorque::init(input_file &inp) {
	BaseForce::init(inp);

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	getInputNumber(&inp, "stiff", &_stiff, 1);
	getInputNumber(&inp, "rate", &_rate, 1);
	getInputNumber(&inp, "base", &_F0, 1);

	std::string strdir;
	double tmpf[3];
	int tmpi;
	getInputString(&inp, "axis", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf, %lf, %lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) {
		throw oxDNAException("Could not parse axis `%s\' for ConstantRateTorque. Aborting", strdir.c_str());
	}
	_axis = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	_axis.normalize();

	getInputString(&inp, "pos0", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf, %lf, %lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) {
		throw oxDNAException("Could not parse pos0 `%s\' for ConstantRateTorque. Aborting", strdir.c_str());
	}
	_pos0 = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);

	getInputString(&inp, "center", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf, %lf, %lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) {
		throw oxDNAException("Could not parse center `%s\' for ConstantRateTorque. Aborting", strdir.c_str());
	}
	_center = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);

	if(getInputString(&inp, "mask", strdir, 0) == KEY_FOUND) {
		tmpi = sscanf(strdir.c_str(), "%lf, %lf, %lf", tmpf, tmpf + 1, tmpf + 2);
		if(tmpi != 3)
			throw oxDNAException("Could not parse mask `%s\' for ConstantRateTorque. Aborting", strdir.c_str());
		_mask = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	}
	else {
		_mask = LR_vector(0., 0., 0.);
	}

	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "ConstantRateTorque");
	std::string description = Utils::sformat("ConstantRateTorque (F0==%g, rate=%g, pos0=%g,%g,%g, axis=%g,%g,%g, center=%g,%g,%g, mask=%g,%g,%g on particle %d", _F0, _rate, _pos0.x, _pos0.y, _pos0.z, _axis.x, _axis.y, _axis.z, _center.x, _center.y, _center.z, _mask.x, _mask.y, _mask.z);

	return std::make_tuple(particle_ids, description);
}

number ConstantRateTorque::potential(llint step, LR_vector &pos) {
	number t = (_F0 + _rate * (number) step);

	number sintheta = sin(t);
	number costheta = cos(t);
	number olcos = ((number) 1.) - costheta;

	number xyo = _axis.x * _axis.y * olcos;
	number xzo = _axis.x * _axis.z * olcos;
	number yzo = _axis.y * _axis.z * olcos;
	number xsin = _axis.x * sintheta;
	number ysin = _axis.y * sintheta;
	number zsin = _axis.z * sintheta;

	LR_matrix R(_axis.x * _axis.x * olcos + costheta, xyo - zsin, xzo + ysin,

	xyo + zsin, _axis.y * _axis.y * olcos + costheta, yzo - xsin,

	xzo - ysin, yzo + xsin, _axis.z * _axis.z * olcos + costheta);

	LR_vector postrap = R * (_pos0 - _center) + _center;

	// we "mask" the resulting vector;
	LR_vector dr = pos - postrap;
	dr.x *= _mask.x;
	dr.y *= _mask.y;
	dr.z *= _mask.z;

	return (number) (0.5 * _stiff * (dr * dr));
}

LR_vector ConstantRateTorque::value(llint step, LR_vector &pos) {
	//
	// dobbiamo ruotare pos0 di (base + rate * t) radianti
	// intorno all'asse passante per il centro.
	// e quello ci da la pos della trappola;
	//
	number t = (_F0 + _rate * (number) step);

	number sintheta = sin(t);
	number costheta = cos(t);
	number olcos = ((number) 1.) - costheta;

	number xyo = _axis.x * _axis.y * olcos;
	number xzo = _axis.x * _axis.z * olcos;
	number yzo = _axis.y * _axis.z * olcos;
	number xsin = _axis.x * sintheta;
	number ysin = _axis.y * sintheta;
	number zsin = _axis.z * sintheta;

	LR_matrix R(_axis.x * _axis.x * olcos + costheta, xyo - zsin, xzo + ysin,

	xyo + zsin, _axis.y * _axis.y * olcos + costheta, yzo - xsin,

	xzo - ysin, yzo + xsin, _axis.z * _axis.z * olcos + costheta);

	LR_vector postrap = R * (_pos0 - _center) + _center;

	// we "mask" the resulting vector;
	number x = -_stiff * (pos.x - postrap.x) * _mask.x;
	number y = -_stiff * (pos.y - postrap.y) * _mask.y;
	number z = -_stiff * (pos.z - postrap.z) * _mask.z;

	return LR_vector(x, y, z);
}
