/*
 * RepulsiveSphere.cpp
 *
 *  Created on: 28/nov/2014
 *      Author: Flavio 
 */

#include "RepulsiveSphere.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

RepulsiveSphere::RepulsiveSphere() :
				BaseForce() {
	_r0 = -1.;
	_r_ext = 1e10;
	_center = LR_vector(0., 0., 0.);
	_rate = 0.;
	_box_ptr = NULL;
}

void RepulsiveSphere::get_settings(input_file &inp) {
	getInputNumber(&inp, "stiff", &this->_stiff, 1);
	getInputNumber(&inp, "r0", &_r0, 1);
	getInputNumber(&inp, "rate", &_rate, 0);
	getInputNumber(&inp, "r_ext", &_r_ext, 0);
	getInputString(&inp, "particle", _particles_string, 1);

	std::string strdir;
	if(getInputString(&inp, "center", strdir, 0) == KEY_FOUND) {
		double tmpf[3];
		int tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
		if(tmpi != 3) throw oxDNAException("Could not parse center %s in external forces file. Aborting", strdir.c_str());
		this->_center = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	}
}

void RepulsiveSphere::init(std::vector<BaseParticle *> & particles, BaseBox *box_ptr) {
	std::string force_description = Utils::sformat("RepulsiveSphere (stiff=%g, r0=%g, rate=%g, center=%g,%g,%g)", this->_stiff, this->_r0, this->_rate, this->_center.x, this->_center.y, this->_center.z);
	this->_add_self_to_particles(particles, _particles_string, force_description);

	_box_ptr = box_ptr;
}

LR_vector RepulsiveSphere::value(llint step, LR_vector &pos) {
	LR_vector dist = _box_ptr->min_image(this->_center, pos);
	number mdist = dist.module();
	number radius = _r0 + _rate * (number) step;

	if(mdist <= radius || mdist >= _r_ext) return LR_vector(0., 0., 0.);
	else return dist * (-this->_stiff * (1. - radius / mdist));
}

number RepulsiveSphere::potential(llint step, LR_vector &pos) {
	LR_vector dist = _box_ptr->min_image(this->_center, pos);
	number mdist = dist.module();
	number radius = _r0 + _rate * (number) step;

	if(mdist <= radius || mdist >= _r_ext) return 0.;
	else return 0.5 * this->_stiff * (mdist - radius) * (mdist - radius);
}
