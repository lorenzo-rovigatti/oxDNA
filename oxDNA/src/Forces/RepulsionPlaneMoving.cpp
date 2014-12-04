/*
 * RepulsionPlaneMoving.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio 
 */

#include "RepulsionPlaneMoving.h"
#include "../Particles/BaseParticle.h"

template<typename number>
RepulsionPlaneMoving<number>::RepulsionPlaneMoving() : BaseForce<number>() {
	_particle = -1;
	_ref_id = -1;
}

template<typename number>
void RepulsionPlaneMoving<number>::get_settings (input_file &inp) {
	getInputInt (&inp, "particle", &_particle, 0);

	getInputNumber(&inp, "stiff", &this->_stiff, 1);
	getInputInt(&inp, "ref_particle", &_ref_id, 1);

	int tmpi;
	double tmpf[3];
	std::string strdir;
	getInputString (&inp, "dir", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if (tmpi != 3) throw oxDNAException ("Could not parse dir %s in external forces file. Aborting", strdir.c_str());
	this->_direction = LR_vector<number> ((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	this->_direction.normalize();
}

template<typename number>
void RepulsionPlaneMoving<number>::init (BaseParticle<number> ** particles, int N, number * my_box_side_ptr) {
	if (_ref_id < 0 || _ref_id >= N) throw oxDNAException ("Trying to add a RepulsionPlaneMoving to non-existent particle %d", _ref_id);
	_p_ptr = particles[_ref_id];

	this->box_side_ptr = my_box_side_ptr;
	
	if (_particle >= N || N < -1) throw oxDNAException ("Trying to add a RepulsionPlaneMoving on non-existent particle %d. Aborting", _particle);
	if (_particle != -1) {
		OX_LOG(Logger::LOG_INFO, "Adding RepulsionPlaneMoving with stiffnes %lf and pos=[%g *x + %g * y +  %g * z + d = 0 ]  on particle %i", this->_stiff,  this->_direction.x, this->_direction.y, this->_direction.z, _particle);
		particles[_particle]->add_ext_force(this);
	}
	else { // force affects all particles
		OX_LOG(Logger::LOG_INFO, "Adding RepulsionPlaneMoving with stiffnes %lf and pos=[%g *x + %g * y +  %g * z + d = 0 ]  on ALL particles", this->_stiff,  this->_direction.x, this->_direction.y, this->_direction.z);
		for (int i = 0; i < N; i ++) particles[i]->add_ext_force(this);
	}
}

template<typename number>
LR_vector<number> RepulsionPlaneMoving<number>::value(llint step, LR_vector<number> &pos) {
	number distance = (pos - _p_ptr->get_abs_pos(*(this->box_side_ptr))) * this->_direction;
	if (distance >=  0)
		return LR_vector<number>(0,0,0);
	else
		return -distance * this->_stiff * this->_direction;
}

template<typename number>
number RepulsionPlaneMoving<number>::potential (llint step, LR_vector<number> &pos) {
	number distance = (pos - _p_ptr->get_abs_pos(*(this->box_side_ptr))) * this->_direction;
//	number distance = _p_ptr->pos * this->_direction; //distance from the plane
	if (distance >= 0)
		return 0;
	else
		return (number) (0.5 * this->_stiff * distance * distance);
}

template class RepulsionPlaneMoving<double>;
template class RepulsionPlaneMoving<float>;
