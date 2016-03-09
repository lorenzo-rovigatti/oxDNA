/*
 * MovingTrap.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio 
 */

#include "MovingTrap.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

template<typename number>
MovingTrap<number>::MovingTrap() : BaseForce<number>() {

}

template<typename number>
void MovingTrap<number>::get_settings (input_file &inp) {
	getInputInt (&inp, "particle", &_particle, 1);

	getInputNumber(&inp, "stiff", &this->_stiff, 1);
	getInputNumber(&inp, "rate", &this->_rate, 1);

	int tmpi;
	double tmpf[3];
	std::string strdir;
	getInputString (&inp, "dir", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if (tmpi != 3) throw oxDNAException ("Could not parse dir %s in external forces file. Aborting", strdir.c_str());
	this->_direction = LR_vector<number> ((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	this->_direction.normalize();

	getInputString (&inp, "pos0", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if (tmpi != 3) throw oxDNAException ("Could not parse pos0 %s in external forces file. Aborting", strdir.c_str());
	this->_pos0 = LR_vector<number> ((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
}

template<typename number>
void MovingTrap<number>::init (BaseParticle<number> ** particles, int N, number * box_side_ptr) {
	if (_particle >= N || N < -1) throw oxDNAException ("Trying to add a MovingTrap on non-existent particle %d. Aborting", _particle);
	if (_particle != -1) {
		OX_LOG (Logger::LOG_INFO, "Adding MovingTrap (stiff=%g, rate=%g, dir=%g,%g,%g, pos0=%g,%g,%g on particle %d", this->_stiff, this->_rate, this->_direction.x, this->_direction.y, this->_direction.z, this->_pos0.x, this->_pos0.y, this->_pos0.z, _particle);
		particles[_particle]->add_ext_force(this);
	}
	else { // force affects all particles
		OX_LOG (Logger::LOG_INFO, "Adding MovingTrap (stiff==%g, rate=%g, dir=%g,%g,%g, pos0=%g,%g,%g on ALL particles", this->_F0, this->_rate, this->_direction.x, this->_direction.y, this->_direction.z, this->_pos0.x, this->_pos0.y, this->_pos0.z);
		for (int i = 0; i < N; i ++) particles[i]->add_ext_force(this);
	}
}

template<typename number>
LR_vector<number> MovingTrap<number>::value(llint step, LR_vector<number> &pos) {
    LR_vector<number> postrap;
    number x, y, z;

    postrap.x = this->_pos0.x + (this->_rate * step) * this-> _direction.x;
    postrap.y = this->_pos0.y + (this->_rate * step) * this-> _direction.y;
    postrap.z = this->_pos0.z + (this->_rate * step) * this-> _direction.z;

    x = - this->_stiff * (pos.x - postrap.x);
    y = - this->_stiff * (pos.y - postrap.y);
    z = - this->_stiff * (pos.z - postrap.z);

    return LR_vector<number>(x, y, z);
}

template<typename number>
number MovingTrap<number>::potential (llint step, LR_vector<number> &pos) {
    LR_vector<number> postrap;

    postrap.x = this->_pos0.x + (this->_rate * step) * this-> _direction.x;
    postrap.y = this->_pos0.y + (this->_rate * step) * this-> _direction.y;
    postrap.z = this->_pos0.z + (this->_rate * step) * this-> _direction.z;

    return (number) (0.5 * this->_stiff * (pos - postrap).norm());
}

template class MovingTrap<double>;
template class MovingTrap<float>;
