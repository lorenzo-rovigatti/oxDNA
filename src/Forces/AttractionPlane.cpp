/**
 * @file    RepulsionPlane.cpp
 * @date    01/aug/2024
 * @author  Matthies, Tilibit 
 *
 */

#include "AttractionPlane.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

AttractionPlane::AttractionPlane() :
				BaseForce() {
	_particle = -1;
	_position = -1.;
}

void AttractionPlane::get_settings (input_file &inp) {
	getInputInt (&inp, "particle", &_particle, 1);

	getInputNumber(&inp, "stiff", &this->_stiff, 1);
	getInputNumber(&inp, "position", &this->_position, 1);

	int tmpi;
	double tmpf[3];
	std::string strdir;
	getInputString (&inp, "dir", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) throw oxDNAException ("Could not parse dir %s in external forces file. Aborting", strdir.c_str());
	this->_direction = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	this->_direction.normalize();
}

void AttractionPlane::init (BaseParticle ** particles, int N, BaseBox * box_ptr) {
	if (_particle >= N || N < -1) throw oxDNAException ("Trying to add a AttractionPlane on non-existent particle %d. Aborting", _particle);
	if (_particle != -1) {
		OX_LOG (Logger::LOG_INFO, "Adding AttractionPlane (stiff=%g, position=%g, dir=%g,%g,%g, on particle %d", this->_stiff, this->_position, this->_direction.x, this->_direction.y, this->_direction.z, _particle);
		particles[_particle]->add_ext_force(this);
	}
	else { // force affects all particles
		OX_LOG (Logger::LOG_INFO, "Adding AttractionPlane (stiff=%g, position=%g, dir=%g,%g,%g, on ALL particles", this->_stiff, this->_position, this->_direction.x, this->_direction.y, this->_direction.z);
		for(int i = 0; i < N; i ++) particles[i]->add_ext_force(this);
	}
}

LR_vector AttractionPlane::value(llint step, LR_vector &pos) {
	number distance = this->_direction*pos + this->_position;

	if(distance >=  0.) 
		return -this->_stiff*this->_direction;
	else
		return -(distance*this->_stiff)*this->_direction;
}

number AttractionPlane::potential(llint step, LR_vector &pos) {
	number distance = this->_direction*pos + this->_position; // distance from the plane
	if(distance >= 0.) 
		return  (number) (this->_stiff*distance);
	else 
		return (number) (0.5*this->_stiff*SQR(distance));
}

