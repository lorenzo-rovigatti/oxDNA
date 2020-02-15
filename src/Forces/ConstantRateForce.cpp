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

ConstantRateForce::ConstantRateForce() :
				BaseForce() {
	_particles_string = "\0";
	dir_as_centre = false;
}

ConstantRateForce::~ConstantRateForce() {

}

void ConstantRateForce::get_settings(input_file &inp) {
	std::string particles_string;
	getInputString(&inp, "particle", _particles_string, 1);

	getInputNumber(&inp, "F0", &this->_F0, 1);
	getInputNumber(&inp, "rate", &this->_rate, 1);
	getInputBool(&inp, "dir_as_centre", &this->dir_as_centre, 0);

	string strdir;
	getInputString(&inp, "dir", strdir, 1);
	vector<string> spl = Utils::split(strdir, ',');
	if(spl.size() != 3) throw oxDNAException("Could not parse 'dir' in external_forces_file. Dying badly");

	this->_direction.x = atof(spl[0].c_str());
	this->_direction.y = atof(spl[1].c_str());
	this->_direction.z = atof(spl[2].c_str());

	if(!dir_as_centre) this->_direction.normalize();
}

void ConstantRateForce::init(std::vector<BaseParticle *> &particles, BaseBox *box_ptr) {
	std::string force_description = Utils::sformat("ConstantRateForce (F=%g, rate=%g, dir=%g,%g,%g)", this->_F0, this->_rate, this->_direction.x, this->_direction.y, this->_direction.z);
	this->_add_self_to_particles(particles, _particles_string, force_description);
}

LR_vector ConstantRateForce::value(llint step, LR_vector &pos) {
	LR_vector dir = this->_direction;
	if(dir_as_centre) {
		dir -= pos;
		dir.normalize();
	}
	LR_vector force = (this->_F0 + this->_rate * step) * dir;
	return force;
}

number ConstantRateForce::potential(llint step, LR_vector &pos) {
	number strength = -(this->_F0 + this->_rate * step);
	if(dir_as_centre) {
		number dist = (this->_direction - pos).module();
		return strength * dist;
	}
	return strength * (pos * this->_direction);
}
