/*
 * LowdimMovingTrap.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio 
 */

#include "LowdimMovingTrap.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

LowdimMovingTrap::LowdimMovingTrap() :
				BaseForce() {
	_particle = -2;
}

void LowdimMovingTrap::get_settings(input_file &inp) {
	getInputInt(&inp, "particle", &_particle, 1);

	getInputNumber(&inp, "stiff", &this->_stiff, 1);
	getInputNumber(&inp, "rate", &this->_rate, 1);

	int tmpi, tmpa[3];
	double tmpf[3];
	std::string strdir;

	getInputString(&inp, "pos0", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf, %lf, %lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) throw oxDNAException("Could not parse pos0 `%s\' for LowdimMivingTrap. Aborting", strdir.c_str());
	this->_pos0 = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);

	getInputString(&inp, "visibility", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%d, %d, %d", tmpa, tmpa + 1, tmpa + 2);
	if(tmpi != 3) throw oxDNAException("Could not parse visibility `%s\' for LowdimMivingTrap. Aborting", strdir.c_str());
	this->_visX = tmpa[0] ? true : false;
	this->_visY = tmpa[1] ? true : false;
	this->_visZ = tmpa[2] ? true : false;

	getInputString(&inp, "dir", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf, %lf, %lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) throw oxDNAException("Could not parse dir `%s\' for LowdimMivingTrap. Aborting", strdir.c_str());
	this->_direction = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	this->_direction.normalize();
}

void LowdimMovingTrap::init(std::vector<BaseParticle *> & particles, int N, BaseBox *box_ptr) {
	if(_particle >= N || N < -1) throw oxDNAException("Trying to add a LowdimMovingTrap on non-existent particle %d. Aborting", _particle);
	if(_particle != -1) {
		OX_LOG(Logger::LOG_INFO, "Adding LowdimMovingTrap (stiff=stiffness %lf and pos=[%g,%g,%g] + (%g * t) [%g,%g,%g] and visX=%i visY=%i visZ=%i on particle %d", this->_stiff, this->_pos0.x, this->_pos0.y, this->_pos0.z, this->_rate, this->_direction.x, this->_direction.y, this->_direction.z, this->_visX, this->_visY, this->_visZ, _particle);
		particles[_particle]->add_ext_force(ForcePtr(this));
	}
	else { // force affects all particles
		OX_LOG (Logger::LOG_INFO, "Adding LowdimMovingTrap (stiff=stiffness %lf and pos=[%g,%g,%g] + (%g * t) [%g,%g,%g] and visX=%i visY=%i visZ=%i on ALL particles", this->_stiff, this->_pos0.x, this->_pos0.y, this->_pos0.z, this->_rate, this->_direction.x, this->_direction.y, this->_direction.z, this->_visX, this->_visY, this->_visZ);
		for (int i = 0; i < N; i ++) {
			particles[i]->add_ext_force(ForcePtr(this));
		}
	}
}

LR_vector LowdimMovingTrap::value(llint step, LR_vector &pos) {
	LR_vector postrap;
	number x, y, z;

	postrap.x = (this->_visX) ? this->_pos0.x + (this->_rate * step) * this->_direction.x : pos.x;
	postrap.y = (this->_visY) ? this->_pos0.y + (this->_rate * step) * this->_direction.y : pos.y;
	postrap.z = (this->_visZ) ? this->_pos0.z + (this->_rate * step) * this->_direction.z : pos.z;

	x = -this->_stiff * (pos.x - postrap.x);
	y = -this->_stiff * (pos.y - postrap.y);
	z = -this->_stiff * (pos.z - postrap.z);

	return LR_vector(x, y, z);
}

number LowdimMovingTrap::potential(llint step, LR_vector &pos) {
	LR_vector postrap;

	postrap.x = (this->_visX) ? this->_pos0.x + (this->_rate * step) * this->_direction.x : pos.x;
	postrap.y = (this->_visY) ? this->_pos0.y + (this->_rate * step) * this->_direction.y : pos.y;
	postrap.z = (this->_visZ) ? this->_pos0.z + (this->_rate * step) * this->_direction.z : pos.z;

	return (number) (0.5 * this->_stiff * (pos - postrap).norm());
}
