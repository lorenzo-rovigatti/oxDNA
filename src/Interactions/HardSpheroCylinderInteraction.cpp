/*
 * HardSpheroCylinderInteraction.cpp
 *
 *  Created on: 29/Oct/2013
 *      Author: Flavio 
 */

#include "HardSpheroCylinderInteraction.h"

template<typename number>
HardSpheroCylinderInteraction<number>::HardSpheroCylinderInteraction() : BaseInteraction<number, HardSpheroCylinderInteraction<number> >() {
	this->_int_map[0] = &HardSpheroCylinderInteraction<number>::_hsc_pot;
}

template<typename number>
HardSpheroCylinderInteraction<number>::~HardSpheroCylinderInteraction() {

}

template<typename number>
void HardSpheroCylinderInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);
	char tmps[512];
	getInputString (&inp, "sim_type", (char *)tmps, 1);
	if (strncmp(tmps, "MC", 512)) throw oxDNAException ("Cannot run HardSpheroCylinder with MD");
	
	// length
	float tmpf;
	getInputFloat (&inp, "length", &tmpf, 1);
	_length = (number) tmpf;
	if (_length < 0.) throw oxDNAException ("Cannot run hard spherocylinders with negative lenght");
	
	// rcut is 2. * (radius + length)
	this->_rcut = (number) 1.001 + _length;
	
	OX_LOG(Logger::LOG_INFO, "Initializing HardSpheroCylinder interaction with length %g [_rcut = %g]", _length, this->_rcut);
}

template<typename number>
void HardSpheroCylinderInteraction<number>::init() {
	this->_sqr_rcut = SQR(this->_rcut);
}

template<typename number>
void HardSpheroCylinderInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new BaseParticle<number>();
}

template<typename number>
number HardSpheroCylinderInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number HardSpheroCylinderInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return (number) 0.f;
}

template<typename number>
number HardSpheroCylinderInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}
	
	return _hsc_pot (p, q, r, update_forces);
}

template<typename number>
void HardSpheroCylinderInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template class HardSpheroCylinderInteraction<float>;
template class HardSpheroCylinderInteraction<double>;
