/*
 * HardWall.cpp
 *
 *  Created on: 5/oct/2016
 *      Author: Lorenzo
 */

#include "HardWall.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

HardWall::HardWall() :
				BaseForce() {
	_particle = -1;
	_position = -1.;
	this->_stiff = 1;
	_sigma = 1.;
}

void HardWall::get_settings(input_file &inp) {
	getInputInt(&inp, "particle", &_particle, 1);

	getInputNumber(&inp, "position", &this->_position, 1);
	getInputNumber(&inp, "sigma", &this->_sigma, 0);

	int tmpi;
	double tmpf[3];
	std::string strdir;
	getInputString(&inp, "dir", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) throw oxDNAException("Could not parse dir %s in external forces file. Aborting", strdir.c_str());
	this->_direction = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	this->_direction.normalize();
}

void HardWall::init(std::vector<BaseParticle *> & particles, BaseBox *box_ptr) {
	int N = particles.size();
	if(_particle >= N || N < -1) throw oxDNAException("Trying to add a HardWall on non-existent particle %d. Aborting", _particle);
	if(_particle != -1) {
		OX_LOG(Logger::LOG_INFO, "Adding HardWall (position=%g, dir=%g,%g,%g, sigma=%g) on particle %d", this->_position, this->_direction.x, this->_direction.y, this->_direction.z, _sigma, _particle);
		particles[_particle]->add_ext_force(ForcePtr(this));
	}
	else { // force affects all particles
		OX_LOG (Logger::LOG_INFO, "Adding HardWall (position=%g, dir=%g,%g,%g, sigma=%g) on ALL particles", this->_position, this->_direction.x, this->_direction.y, this->_direction.z, _sigma);
		for (int i = 0; i < N; i ++) {
			particles[i]->add_ext_force(ForcePtr(this));
		}
	}
}

LR_vector HardWall::value(llint step, LR_vector &pos) {
	throw oxDNAException("HardWall can be used only in Monte Carlo simulations.");
}

number HardWall::potential(llint step, LR_vector &pos) {
	number distance = (this->_direction * pos + this->_position) / _sigma; // distance from the plane in units of _sigma
	if(distance < 1.) return 10e8;
	return 0.;
}
