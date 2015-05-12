/*
 * HardCylinderInteraction.cpp
 *
 *  Created on: 31/Oct/2013
 *      Author: Flavio 
 */

#include "HardCylinderInteraction.h"

template<typename number>
HardCylinderInteraction<number>::HardCylinderInteraction() : BaseInteraction<number, HardCylinderInteraction<number> >() {
	this->_int_map[HardCylinder] = &HardCylinderInteraction<number>::_hc_pot;
}

template<typename number>
HardCylinderInteraction<number>::~HardCylinderInteraction() {

}

template<typename number>
void HardCylinderInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);
	char tmps[512];
	getInputString (&inp, "sim_type", (char *)tmps, 1);
	if (strncmp(tmps, "MC", 512)) throw oxDNAException ("Cannot run Hard Cylinders with MD");
	
	float tmpf;
	getInputFloat (&inp, "height", &tmpf, 1);
	_height = (number) tmpf;
	
	OX_LOG(Logger::LOG_INFO, "Initializing HardCylinder interaction with height %g", _height);
	
	// r_cut for outer bounding box
	this->_rcut = (number) 1.001 * 2. * sqrt (0.5 * 0.5  + (_height / 2.) * (_height / 2.));

	OX_LOG(Logger::LOG_INFO, "Using r_cut of %g", this->_rcut);
}

template<typename number>
void HardCylinderInteraction<number>::init() {
	this->_sqr_rcut = SQR(this->_rcut);
}

template<typename number>
void HardCylinderInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new BaseParticle<number>();
}

template<typename number>
number HardCylinderInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number HardCylinderInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return (number) 0.f;
}

template<typename number>
number HardCylinderInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}
	
	return _hc_pot (p, q, r, update_forces);
}

template<typename number>
void HardCylinderInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template class HardCylinderInteraction<float>;
template class HardCylinderInteraction<double>;

