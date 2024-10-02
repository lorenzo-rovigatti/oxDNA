/*
 * AlignmentField.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio 
 */

#include "AlignmentField.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

AlignmentField::AlignmentField() :
				BaseForce() {
	_particle = -1;
	_F = -1.f;
	_v_idx = -1;
	_v_ptr = NULL;
}

std::tuple<std::vector<int>, std::string> AlignmentField::init(input_file &inp) {
	BaseForce::init(inp);

	getInputInt(&inp, "particle", &_particle, 0);
	getInputInt(&inp, "v_idx", &_v_idx, 1);
	getInputNumber(&inp, "F", &_F, 1);
	if(_v_idx < 0 || _v_idx >= 6) throw oxDNAException("(AlignmentField.cpp) v_idx must be >= 0 and <= 5, got %d. Aborting", _v_idx);

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

	std::vector<BaseParticle *> &particles = CONFIG_INFO->particles();

	int N = particles.size();
	if(_particle >= N || N < 0) throw oxDNAException("Trying to add a AlignmentField on non-existent particle %d. Aborting", _particle);
	switch(_v_idx) {
	case 0:
		_v_ptr = &(particles[_particle]->orientation.v1);
		break;
	case 1:
		_v_ptr = &(particles[_particle]->orientation.v2);
		break;
	case 2:
		_v_ptr = &(particles[_particle]->orientation.v3);
		break;
	case 3:
		_v_ptr = &(particles[_particle]->orientationT.v1);
		break;
	case 4:
		_v_ptr = &(particles[_particle]->orientationT.v2);
		break;
	case 5:
		_v_ptr = &(particles[_particle]->orientationT.v3);
		break;
	default:
		throw oxDNAException("Should Never Get here %s %s", __FILE__, __LINE__);
		break;
	}

	std::string description = Utils::sformat("AlignmentField (F=%g, dir=%g,%g,%g)", _F, _direction.x, _direction.y, _direction.z);
	return std::make_tuple(std::vector<int>{_particle}, description);
}

LR_vector AlignmentField::value(llint step, LR_vector &pos) {
	throw oxDNAException("Not implemented... %s %s", __FILE__, __LINE__);
}

number AlignmentField::potential(llint step, LR_vector &pos) {
	return _F * (1.f - ((*_v_ptr) * (_direction)));
}
