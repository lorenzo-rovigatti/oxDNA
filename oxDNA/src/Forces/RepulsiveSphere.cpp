/*
 * RepulsiveSphere.cpp
 *
 *  Created on: 28/nov/2014
 *      Author: Flavio 
 */

#include "RepulsiveSphere.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

template<typename number>
RepulsiveSphere<number>::RepulsiveSphere() : BaseForce<number>() {
	_particle = -1;
	_r0 = -1.;
	_center = LR_vector<number>(0., 0., 0.);
	_rate = 0.;
	_box_side_ptr = NULL;
}

template<typename number>
void RepulsiveSphere<number>::get_settings (input_file &inp) {
	getInputNumber(&inp, "stiff", &this->_stiff, 1);
	getInputNumber(&inp, "r0", &_r0, 1);
	getInputNumber(&inp, "rate", &_rate, 0);

	std::string strdir;
	if (getInputString (&inp, "center", strdir, 0) == KEY_FOUND) {
		double tmpf[3];
		int tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
		if (tmpi != 3) throw oxDNAException ("Could not parse center %s in external forces file. Aborting", strdir.c_str());
		this->_center = LR_vector<number> ((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	}
}

template<typename number>
void RepulsiveSphere<number>::init (BaseParticle<number> ** particles, int N, number * my_box_side_ptr) {
	if (this->_particle >= N || N < -1) throw oxDNAException ("Trying to add a RepulsiveSphere on non-existent particle %d. Aborting", this->_particle);
	if (this->_particle != -1) {
		OX_LOG (Logger::LOG_INFO, "Adding RepulsiveSphere force (stiff=%g, r0=%g, rate=%g, center=%g,%g,%g) on particle %d", this->_stiff, this->_r0, this->_rate, this->_center.x, this->_center.y, this->_center.z, _particle);
		particles[_particle]->add_ext_force(this);
	}
	else { // force affects all particles
		OX_LOG (Logger::LOG_INFO, "Adding RepulsiveSphere force (stiff=%g, r0=%g, rate=%g, center=%g,%g,%g) on ALL particles", this->_stiff, this->_r0, this->_rate, this->_center.x, this->_center.y, this->_center.z);
		for (int i = 0; i < N; i ++) particles[i]->add_ext_force(this);
	}
	_box_side_ptr = my_box_side_ptr; 
}

template<typename number>
LR_vector<number> RepulsiveSphere<number>::value(llint step, LR_vector<number> &pos) {
	LR_vector<number> dist = pos.minimum_image(this->_center, *(this->_box_side_ptr));
	number mdist = dist.module();
	number radius = _r0 + _rate * (number) step;

	if(mdist <= radius) return LR_vector<number>(0., 0., 0.);
	else return dist * (- this->_stiff * (1. - radius / mdist));
}

template<typename number>
number RepulsiveSphere<number>::potential (llint step, LR_vector<number> &pos) {
	LR_vector<number> dist = pos.minimum_image(this->_center, *(this->_box_side_ptr));
	number mdist = dist.module();
	number radius = _r0 + _rate * (number) step;

	if (mdist <= radius) return 0.;
	else return 0.5 * this->_stiff * (mdist - radius) * (mdist - radius);
}

template class RepulsiveSphere<double>;
template class RepulsiveSphere<float>;

