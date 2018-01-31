/*
 * ConstantRateForce.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio 
 *              Modified by Ferdinando on 11/Dec/2015
 */

#include "ConstantRateForce.h"
#include "../Particles/BaseParticle.h"

using namespace std;

template<typename number>
ConstantRateForce<number>::ConstantRateForce() : BaseForce<number>() {
	_particles_string = "\0";
	dir_as_centre = false;
}

template<typename number>
ConstantRateForce<number>::~ConstantRateForce() {

}

template<typename number>
void ConstantRateForce<number>::get_settings(input_file &inp) {
	std::string particles_string;
	getInputString (&inp, "particle", _particles_string, 1);

	getInputNumber (&inp, "F0", &this->_F0, 1);
	getInputNumber (&inp, "rate", &this->_rate, 1);
	getInputBool(&inp, "dir_as_centre", &this->dir_as_centre, 0);

	string strdir;
	getInputString (&inp, "dir", strdir, 1);
	vector<string> spl = Utils::split(strdir, ',');
	if(spl.size() != 3) throw oxDNAException("Could not parse 'dir' in external_forces_file. Dying badly");

	this->_direction.x = atof(spl[0].c_str());
	this->_direction.y = atof(spl[1].c_str());
	this->_direction.z = atof(spl[2].c_str());

	if(!dir_as_centre) this->_direction.normalize();
}

template<typename number>
void ConstantRateForce<number>::init(BaseParticle<number> **particles, int N, BaseBox<number> * box_ptr) {
	std::string force_description = Utils::sformat("ConstantRateForce (F=%g, rate=%g, dir=%g,%g,%g)", this->_F0, this->_rate, this->_direction.x, this->_direction.y, this->_direction.z);
	this->_add_self_to_particles(particles, N, _particles_string, force_description);
}

template<typename number>
LR_vector<number> ConstantRateForce<number>::value(llint step, LR_vector<number> &pos) {
	LR_vector<number> dir = this->_direction;
	if(dir_as_centre) {
		dir -= pos;
		dir.normalize();
	}
	LR_vector<number> force = (this->_F0 + this->_rate*step)*dir;
	return force;
}

template<typename number>
number ConstantRateForce<number>::potential(llint step, LR_vector<number> &pos) {
	number strength = -(this->_F0 + this->_rate * step);
	if(dir_as_centre) {
		number dist = (this->_direction - pos).module();
		return strength * dist;
	}
	return strength * (pos * this->_direction);
}

template class ConstantRateForce<double>;
template class ConstantRateForce<float>;

