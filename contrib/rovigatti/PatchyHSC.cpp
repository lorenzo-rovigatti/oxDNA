/*
 * PatchyHSC.cpp
 *
 *  Created on: 11/ago/2014
 *      Author: lorenzo
 */

#include "PatchyHSC.h"

template <typename number>
PatchyHSC<number>::PatchyHSC() : HardSpheroCylinderInteraction<number>() {
	//this->_int_map[PATCHY] = &PatchyHSC<number>::_patchy;
}

template <typename number>
PatchyHSC<number>::~PatchyHSC() {

}

template<typename number>
void PatchyHSC<number>::get_settings(input_file &inp) {
	HardSpheroCylinderInteraction<number>::get_settings(inp);

	getInputNumber(&inp, "PHSC_protrusion", &_protrusion, 1);
	getInputNumber(&inp, "PHSC_patch_r", &_patch_r, 1);
}

template<typename number>
void PatchyHSC<number>::init() {
	HardSpheroCylinderInteraction<number>::init();

	_centre_patch_dist = 0.5*this->_rcut + _protrusion;
	_sqr_patch_rcut = SQR(2*_patch_r);

	this->_rcut = this->_rcut + 2*_protrusion + 2*_patch_r;
	this->_sqr_rcut = SQR(this->_rcut);

	OX_LOG(Logger::LOG_INFO, "Initialized PatchyHSC interaction with rcut %g", this->_rcut);
}

template<typename number>
number PatchyHSC<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number energy = HardSpheroCylinderInteraction<number>::pair_interaction_nonbonded(p, q, r, update_forces);
	energy += _patchy(p, q, r, update_forces);

	return energy;
}

template class PatchyHSC<float>;
template class PatchyHSC<double>;
