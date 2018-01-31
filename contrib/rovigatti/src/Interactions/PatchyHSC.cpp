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
void PatchyHSC<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new PatchySpherocylinder<number>(_centre_patch_dist);
}

template<typename number>
void PatchyHSC<number>::init() {
	HardSpheroCylinderInteraction<number>::init();

	_spherocylinder_length = (number) 1.001 + this->_length;
	_sqr_spherocylinder_length = SQR(_spherocylinder_length);

	_centre_patch_dist = 0.5*_spherocylinder_length + _protrusion;
	_sqr_patch_rcut = SQR(2*_patch_r);
	_sqr_patch_shoulder_rcut = SQR(0.5*_patch_r);

	this->_rcut = this->_rcut + 2*_protrusion + 2*_patch_r;
	this->_sqr_rcut = SQR(this->_rcut);

	OX_LOG(Logger::LOG_INFO, "Initialized PatchyHSC interaction with rcut %g", this->_rcut);
}

template<typename number>
number PatchyHSC<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(this->_box_side > 0.1 && this->_box_side < 2*this->_rcut) throw oxDNAException("The box should be larger than twice the effective diameter of the particles (%lf)\n", 2*this->_rcut);

	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	number sqr_r = r->norm();

	if(sqr_r < _sqr_spherocylinder_length && InteractionUtils::spherocylinder_overlap(*r, p->orientation.v3, q->orientation.v3, this->_length)) {
		this->set_is_infinite(true);
		return 1.0e12;
	}
	if(sqr_r < this->_sqr_rcut) return _patchy(p, q, r, update_forces);

	return 0.;
}

template class PatchyHSC<float>;
template class PatchyHSC<double>;
