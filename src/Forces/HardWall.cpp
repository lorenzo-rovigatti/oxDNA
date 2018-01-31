/*
 * HardWall.cpp
 *
 *  Created on: 5/oct/2016
 *      Author: Lorenzo
 */

#include "HardWall.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

template<typename number>
HardWall<number>::HardWall() : BaseForce<number>() {
	_particle = -1;
	_position = -1.;
	this->_stiff = 1;
	_sigma = 1.;
}

template<typename number>
void HardWall<number>::get_settings (input_file &inp) {
	getInputInt (&inp, "particle", &_particle, 1);

	getInputNumber(&inp, "position", &this->_position, 1);
	getInputNumber(&inp, "sigma", &this->_sigma, 0);

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
void HardWall<number>::init (BaseParticle<number> ** particles, int N, BaseBox<number> * box_ptr) {
	if (_particle >= N || N < -1) throw oxDNAException ("Trying to add a HardWall on non-existent particle %d. Aborting", _particle);
	if (_particle != -1) {
		OX_LOG (Logger::LOG_INFO, "Adding HardWall (position=%g, dir=%g,%g,%g, sigma=%g) on particle %d", this->_position, this->_direction.x, this->_direction.y, this->_direction.z, _sigma, _particle);
		particles[_particle]->add_ext_force(this);
	}
	else { // force affects all particles
		OX_LOG (Logger::LOG_INFO, "Adding HardWall (position=%g, dir=%g,%g,%g, sigma=%g) on ALL particles", this->_position, this->_direction.x, this->_direction.y, this->_direction.z, _sigma);
		for (int i = 0; i < N; i ++) particles[i]->add_ext_force(this);
	}
}

template<typename number>
LR_vector<number> HardWall<number>::value(llint step, LR_vector<number> &pos) {
	throw oxDNAException("HardWall can be used only in Monte Carlo simulations.");
}

template<typename number>
number HardWall<number>::potential (llint step, LR_vector<number> &pos) {
	number distance = (this->_direction*pos + this->_position)/_sigma; // distance from the plane in units of _sigma
	if(distance < 1.) return 10e8;
	return 0.;
}

template class HardWall<double>;
template class HardWall<float>;
