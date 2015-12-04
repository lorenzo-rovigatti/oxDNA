/*
 * RepulsionPlane.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio 
 */

#include "RepulsionPlane.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

template<typename number>
RepulsionPlane<number>::RepulsionPlane() : BaseForce<number>() {
	_particle = -1;
	_position = -1.;
}

template<typename number>
void RepulsionPlane<number>::get_settings (input_file &inp) {
	getInputInt (&inp, "particle", &_particle, 0);

	getInputNumber(&inp, "stiff", &this->_stiff, 1);
	getInputNumber(&inp, "position", &this->_position, 1);

	int tmpi;
	double tmpf[3];
	std::string strdir;
	getInputString (&inp, "dir", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) throw oxDNAException ("Could not parse dir %s in external forces file. Aborting", strdir.c_str());
	this->_direction = LR_vector<number> ((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	this->_direction.normalize();
}

template<typename number>
void RepulsionPlane<number>::init (BaseParticle<number> ** particles, int N, number * box_side_ptr) {
	if (_particle >= N || N < -1) throw oxDNAException ("Trying to add a RepulsionPlane on non-existent particle %d. Aborting", _particle);
	if (_particle != -1) {
		OX_LOG (Logger::LOG_INFO, "Adding RepulsionPlane (stiff=%g, position=%g, dir=%g,%g,%g, on particle %d", this->_stiff, this->_position, this->_direction.x, this->_direction.y, this->_direction.z, _particle);
		particles[_particle]->add_ext_force(this);
	}
	else { // force affects all particles
		OX_LOG (Logger::LOG_INFO, "Adding RepulsionPlane (stiff=%g, position=%g, dir=%g,%g,%g, on ALL particles", this->_stiff, this->_position, this->_direction.x, this->_direction.y, this->_direction.z);
		for(int i = 0; i < N; i ++) particles[i]->add_ext_force(this);
	}
}

template<typename number>
LR_vector<number> RepulsionPlane<number>::value(llint step, LR_vector<number> &pos) {
	number distance = this->_direction*pos + this->_position;
	if(distance >=  0.) return LR_vector<number>(0., 0., 0.);
	else return -(distance*this->_stiff)*this->_direction;
}

template<typename number>
number RepulsionPlane<number>::potential(llint step, LR_vector<number> &pos) {
	number distance = this->_direction*pos + this->_position; // distance from the plane
	if(distance >= 0.) return 0.;
	else return (number) (0.5*this->_stiff*SQR(distance));
}

template class RepulsionPlane<double>;
template class RepulsionPlane<float>;
