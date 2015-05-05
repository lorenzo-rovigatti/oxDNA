/*
 * ExternalForce.cpp
 *
 *  Created on: 8/apr/2015
 *      Author: Megan
 */

#include "SawtoothForce.h"
#include "../Particles/BaseParticle.h"


template<typename number>
SawtoothForce<number>::SawtoothForce() : BaseForce<number>() {
	_particle = 0;
	_wait_time = 0.0;
	_increment = 0.0;
}

template<typename number>
SawtoothForce<number>::~SawtoothForce() {

}

template<typename number>
void SawtoothForce<number>::get_settings(input_file &inp) {
	getInputInt (&inp, "particle", &this->_particle, 1);

	getInputNumber (&inp, "F0", &this->_F0, 1);
	getInputNumber (&inp, "wait_time", &this->_wait_time, 1);
	getInputNumber (&inp, "increment", &this->_increment, 1);

	std::string strdir;
	getInputString (&inp, "dir", strdir, 1);
	double x, y, z;
	int tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", &x, &y, &z);
	if (tmpi != 3) throw oxDNAException("could not parse direction in external_forces_file. Dying badly");
	this->_direction = LR_vector<number> ((number) x, (number) y, number (z));
	this->_direction.normalize();
}

template<typename number>
void SawtoothForce<number>::init(BaseParticle<number> **particles, int N, number *box_side) {
	if (_particle >= N || N < -1) throw oxDNAException ("Trying to add a SawtoothForce on non-existent particle %d. Aborting", _particle);
	if (_particle != -1) {
		OX_LOG (Logger::LOG_INFO, "Adding SawtoothForce (F==%g, wait_time=%g, increment=%g, dir=%g,%g,%g on particle %d", this->_F0, this->_wait_time, this->_increment, this->_direction.x, this->_direction.y, this->_direction.z, _particle);
		particles[_particle]->add_ext_force(this);
	}
	else { // force affects all particles
		OX_LOG (Logger::LOG_INFO, "Adding SawtoothForce (F==%g, wait_time=%g, increment=%g, dir=%g,%g,%g on ALL particles", this->_F0, this->_wait_time, this->_increment, this->_direction.x, this->_direction.y, this->_direction.z);
		for (int i = 0; i < N; i ++) particles[i]->add_ext_force(this);
	}
}

template<typename number>
LR_vector<number> SawtoothForce<number>::value(llint step, LR_vector<number> &pos) {
    number x = (this->_F0  + (int)((step-1)/this->_wait_time) * this->_increment) * this->_direction.x;
    number y = (this->_F0  + (int)((step-1)/this->_wait_time) * this->_increment) * this->_direction.y;
    number z = (this->_F0  + (int)((step-1)/this->_wait_time) * this->_increment) * this->_direction.z;
    return LR_vector<number>(x, y, z);
}

template<typename number>
number SawtoothForce<number>::potential(llint step, LR_vector<number> &pos) {
  return (number) -(this->_F0  + (int)((step-1)/this->_wait_time) * this->_increment) * (pos * this->_direction);
}

template class SawtoothForce<double>;
template class SawtoothForce<float>;

