/*
 * RepulsiveSphere.cpp
 *
 *  Created on: 28/nov/2014
 *      Author: Flavio 
 */

#include "RepulsiveSphereSmooth.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

template<typename number>
RepulsiveSphereSmooth<number>::RepulsiveSphereSmooth() :
				BaseForce<number>() {
	_particle = -1;
	_r0 = -1.;
	_r_ext = 1e10;
	_center = LR_vector<number>(0., 0., 0.);
	_alpha = -1.;
	_smooth = -1.;
	_box_ptr = NULL;
}

template<typename number>
void RepulsiveSphereSmooth<number>::get_settings(input_file &inp) {
	getInputNumber(&inp, "stiff", &this->_stiff, 1);
	getInputNumber(&inp, "r0", &_r0, 1);
	getInputNumber(&inp, "r_ext", &_r_ext, 0);
	getInputNumber(&inp, "r_ext", &_smooth, 1);
	getInputNumber(&inp, "r_ext", &_alpha, 1);
	getInputInt(&inp, "particle", &_particle, 1);

	std::string strdir;
	if(getInputString(&inp, "center", strdir, 0) == KEY_FOUND) {
		double tmpf[3];
		int tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
		if(tmpi != 3) throw oxDNAException("Could not parse center %s in external forces file. Aborting", strdir.c_str());
		this->_center = LR_vector<number>((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	}
}

template<typename number>
void RepulsiveSphereSmooth<number>::init(BaseParticle<number> ** particles, int N, BaseBox<number> * box_ptr) {
	if(this->_particle >= N || N < -1) throw oxDNAException("Trying to add a RepulsiveSphereSmooth on non-existent particle %d. Aborting", this->_particle);
	if(this->_particle != -1) {
		OX_LOG(Logger::LOG_INFO, "Adding RepulsiveSphereSmooth force (stiff=%g, r0=%g,  center=%g,%g,%g) on particle %d", this->_stiff, this->_r0, this->_center.x, this->_center.y, this->_center.z, _particle);
		particles[_particle]->add_ext_force(this);
	}
	else { // force affects all particles
		OX_LOG (Logger::LOG_INFO, "Adding RepulsiveSphereSmooth force (stiff=%g, r0=%g,  center=%g,%g,%g) on ALL particles", this->_stiff, this->_r0, this->_center.x, this->_center.y, this->_center.z);
		for (int i = 0; i < N; i ++) particles[i]->add_ext_force(this);
	}
	_box_ptr = box_ptr;
}

template<typename number>
LR_vector<number> RepulsiveSphereSmooth<number>::value(llint step, LR_vector<number> &pos) {
	LR_vector<number> dist = _box_ptr->min_image(this->_center, pos);
	number mdist = dist.module();

	if(mdist < _r0 || mdist > _r_ext) {
		return LR_vector<number>(0., 0., 0.);
	}
	else {
		if(mdist >= _alpha && mdist <= _r_ext) {
			return dist * (-(this->_stiff * 0.5 * exp((mdist - _alpha) / _smooth)) / mdist);
		}
		else {
			return dist * (-(this->_stiff * mdist - this->_stiff * 0.5 * exp(-(mdist - _alpha) / _smooth)) / mdist);
		}
	}
}

template<typename number>
number RepulsiveSphereSmooth<number>::potential(llint step, LR_vector<number> &pos) {
	LR_vector<number> dist = _box_ptr->min_image(this->_center, pos);
	number mdist = dist.module();

	if(mdist < _r0 || mdist > _r_ext) return 0.;
	else {
		if(mdist >= _alpha && mdist <= _r_ext) return (this->_stiff * 0.5 * _smooth * exp((mdist - _alpha) / _smooth));
		else return (this->_stiff * mdist + this->_stiff * 0.5 * _smooth * exp(-(mdist - _alpha) / _smooth));
	}
}

template class RepulsiveSphereSmooth<double> ;
template class RepulsiveSphereSmooth<float> ;

