/*
 * LJCone.cpp
 *
 *  Created on: 24/aug/2017
 *      Author: Lorenzo
 */

#include "LJCone.h"

#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

template<typename number>
LJCone<number>::LJCone() : BaseForce<number>() {
	_particle = -1;
	this->_stiff = 1;
	_n = 6;
	_sigma = 1.;
	_cutoff = 1e6;
	_generate_inside = false;
	_only_repulsive = false;
	_box = NULL;
}

template<typename number>
void LJCone<number>::get_settings (input_file &inp) {
	getInputInt (&inp, "particle", &_particle, 1);

	getInputNumber(&inp, "stiff", &this->_stiff, 0);
	getInputNumber(&inp, "sigma", &_sigma, 0);
	getInputNumber(&inp, "alpha", &_alpha, 1);
	getInputInt(&inp, "n", &_n, 0);
	if(_n % 2 != 0) throw oxDNAException("RepulsiveCone: n (%d) should be an even integer. Aborting", _n);

	getInputBool(&inp, "only_repulsive", &_only_repulsive, 0);
	if(_only_repulsive) _cutoff = pow(2., 1./_n);

	getInputBool(&inp, "generate_inside", &_generate_inside, 0);

	int tmpi;
	double tmpf[3];
	std::string str_vector;
	getInputString(&inp, "dir", str_vector, 1);
	tmpi = sscanf(str_vector.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if (tmpi != 3) throw oxDNAException ("RepulsiveCone: could not parse dir %s in external forces file. Aborting", str_vector.c_str());
	this->_direction = LR_vector<number> ((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	this->_direction.normalize();

	getInputString(&inp, "pos0", str_vector, 1);
	tmpi = sscanf(str_vector.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) throw oxDNAException ("RepulsiveCone: could not parse pos0 %s in external forces file. Aborting", str_vector.c_str());
	this->_pos0 = LR_vector<number> ((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
}

template<typename number>
void LJCone<number>::init (BaseParticle<number> ** particles, int N, BaseBox<number> *box_ptr) {
	_box = box_ptr;
	_sin_alpha = sin(_alpha);
	_cos_alpha = cos(_alpha);
	_tan_alpha = tan(_alpha);

	if (_particle >= N || N < -1) throw oxDNAException ("Trying to add a RepulsiveCone on non-existent particle %d. Aborting", _particle);
	if (_particle != -1) {
		OX_LOG(Logger::LOG_INFO, "Adding RepulsiveCone (stiff=%g, pos0=%g,%g,%g, dir=%g,%g,%g, sigma=%g, n=%d) on particle %d", this->_stiff, this->_pos0.x, this->_pos0.y, this->_pos0.z, this->_direction.x, this->_direction.y, this->_direction.z, _sigma, _n, _particle);
		particles[_particle]->add_ext_force(this);
	}
	else { // force affects all particles
		OX_LOG(Logger::LOG_INFO, "Adding RepulsiveCone (stiff=%g, pos0=%g,%g,%g, dir=%g,%g,%g, sigma=%g, n=%d) on ALL particles", this->_stiff, this->_pos0.x, this->_pos0.y, this->_pos0.z, this->_direction.x, this->_direction.y, this->_direction.z, _sigma, _n);
		for (int i = 0; i < N; i ++) particles[i]->add_ext_force(this);
	}
}

template<typename number>
LR_vector<number> LJCone<number>::value(llint step, LR_vector<number> &pos) {
	LR_vector<number> v_from_apex = pos - this->_pos0;

	number d_along_axis = v_from_apex*this->_direction;
	LR_vector<number> v_along_axis = this->_direction*d_along_axis;

	LR_vector<number> v_from_axis = v_along_axis - v_from_apex;
	number d_from_axis = v_from_axis.module();
	number d_from_cone = d_along_axis*_sin_alpha - d_from_axis*_cos_alpha;
	number rel_distance = d_from_cone/_sigma;

	if(rel_distance > _cutoff) return LR_vector<number>(0., 0., 0.);

	// now we compute the normal to the cone
	number C = d_from_axis*_tan_alpha;
	LR_vector<number> C_v = (d_along_axis + C)*this->_direction;
	LR_vector<number> normal = C_v - v_from_apex;
	normal.normalize();

	number lj_part = pow(rel_distance, -_n);
	return normal*(4*_n*this->_stiff*(2*SQR(lj_part) - lj_part)/d_from_cone);
}

template<typename number>
number LJCone<number>::potential(llint step, LR_vector<number> &pos) {
	LR_vector<number> v_from_apex = pos - this->_pos0;

	number d_along_axis = v_from_apex*this->_direction;
	LR_vector<number> v_along_axis = this->_direction*d_along_axis;
	LR_vector<number> v_from_axis = v_from_apex - v_along_axis;
	number d_from_cone = d_along_axis*_sin_alpha - v_from_axis.module()*_cos_alpha;
	number rel_distance = d_from_cone/_sigma; // distance from the plane in units of _sigma

	if(_generate_inside && rel_distance < 0.) return 10e8;
	if(rel_distance > _cutoff) rel_distance = _cutoff;
	number lj_part = pow(rel_distance, -_n);
	return 4*this->_stiff*(SQR(lj_part) - lj_part);
}

template class LJCone<double>;
template class LJCone<float>;
