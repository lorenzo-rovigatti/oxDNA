/*
 * LJWall.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio 
 */

#include "LJWall.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

template<typename number>
LJWall<number>::LJWall() : BaseForce<number>() {
	_particle = -1;
	_position = -1.;
	this->_stiff = 1;
	_n = 6;
	_sigma = 1.;
	_cutoff = 1e6;
	_generate_inside = false;
	_only_repulsive = false;
}

template<typename number>
void LJWall<number>::get_settings (input_file &inp) {
	getInputInt (&inp, "particle", &_particle, 0);

	getInputNumber(&inp, "stiff", &this->_stiff, 0);
	getInputNumber(&inp, "position", &this->_position, 1);
	getInputNumber(&inp, "sigma", &this->_sigma, 0);
	getInputInt(&inp, "n", &_n, 0);
	if(_n % 2 != 0) throw oxDNAException("LJWall: n (%d) should be an even integer. Aborting", _n);

	getInputBool(&inp, "only_repulsive", &_only_repulsive, 0);
	if(_only_repulsive) _cutoff = pow(2., 1./_n);

	getInputBool(&inp, "generate_inside", &_generate_inside, 0);

	int tmpi;
	double tmpf[3];
	std::string strdir;
	getInputString (&inp, "dir", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if (tmpi != 3) throw oxDNAException ("Could not parse dir %s in external forces file. Aborting", strdir.c_str());
	this->_direction = LR_vector<number> ((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	this->_direction.normalize();
}

template<typename number>
void LJWall<number>::init (BaseParticle<number> ** particles, int N, number * box_side_ptr) {
	if (_particle >= N || N < -1) throw oxDNAException ("Trying to add a LJWall on non-existent particle %d. Aborting", _particle);
	if (_particle != -1) {
		OX_LOG (Logger::LOG_INFO, "Adding LJWall (stiff=%g, position=%g, dir=%g,%g,%g, sigma=%g, n=%d) on particle %d", this->_stiff, this->_position, this->_direction.x, this->_direction.y, this->_direction.z, _sigma, _n, _particle);
		particles[_particle]->add_ext_force(this);
	}
	else { // force affects all particles
		OX_LOG (Logger::LOG_INFO, "Adding LJWall (stiff=%g, position=%g, dir=%g,%g,%g, sigma=%g, n=%d) on ALL particles", this->_stiff, this->_position, this->_direction.x, this->_direction.y, this->_direction.z, _sigma, _n);
		for (int i = 0; i < N; i ++) particles[i]->add_ext_force(this);
	}
}

template<typename number>
LR_vector<number> LJWall<number>::value(llint step, LR_vector<number> &pos) {
	number distance = this->_direction*pos + this->_position; // distance from the plane
	number rel_distance = distance/_sigma; // distance from the plane in units of _sigma
	if(rel_distance > _cutoff) return LR_vector<number>(0., 0., 0.);
	number lj_part = pow(rel_distance, -_n);
	return this->_direction*(4*_n*this->_stiff*(2*SQR(lj_part) - lj_part)/distance);
}

template<typename number>
number LJWall<number>::potential (llint step, LR_vector<number> &pos) {
	number distance = (this->_direction*pos + this->_position)/_sigma; // distance from the plane in units of _sigma
	if(_generate_inside && distance < 0.) return 10e8;
	if(distance > _cutoff) distance = _cutoff;
	number lj_part = pow(distance, -_n);
	return 4*this->_stiff*(SQR(lj_part) - lj_part);
}

template class LJWall<double>;
template class LJWall<float>;
