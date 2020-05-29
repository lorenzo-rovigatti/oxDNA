/*
 * PatchyHSC.cpp
 *
 *  Created on: 11/ago/2014
 *      Author: lorenzo
 */

#include "PatchyHSC.h"

PatchyHSC::PatchyHSC() :
				HardSpheroCylinderInteraction() {

}

PatchyHSC::~PatchyHSC() {

}

void PatchyHSC::get_settings(input_file &inp) {
	HardSpheroCylinderInteraction::get_settings(inp);

	getInputNumber(&inp, "PHSC_protrusion", &_protrusion, 1);
	getInputNumber(&inp, "PHSC_patch_r", &_patch_r, 1);

	CHECK_BOX("PatchyHSC", inp);
}

void PatchyHSC::allocate_particles(std::vector<BaseParticle *> &particles) {
	for(int i = 0; i < CONFIG_INFO->N(); i++)
		particles[i] = new PatchySpherocylinder(_centre_patch_dist);
}

void PatchyHSC::init() {
	HardSpheroCylinderInteraction::init();

	_spherocylinder_length = (number) 1.001 + _length;
	_sqr_spherocylinder_length = SQR(_spherocylinder_length);

	_centre_patch_dist = 0.5 * _spherocylinder_length + _protrusion;
	_sqr_patch_rcut = SQR(2 * _patch_r);
	_sqr_patch_shoulder_rcut = SQR(0.5 * _patch_r);

	_rcut = _rcut + 2 * _protrusion + 2 * _patch_r;
	_sqr_rcut = SQR(_rcut);

	OX_LOG(Logger::LOG_INFO, "Initialized PatchyHSC interaction with rcut %g", _rcut);
}

number PatchyHSC::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_box->box_sides()[0] > 0.1 && _box->box_sides()[0] < 2 * _rcut) {
		throw oxDNAException("The box should be larger than twice the effective diameter of the particles (%lf)\n", 2 * _rcut);
	}

	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	number sqr_r = _computed_r.norm();

	if(sqr_r < _sqr_spherocylinder_length && InteractionUtils::spherocylinder_overlap(_computed_r, p->orientation.v3, q->orientation.v3, _length)) {
		set_is_infinite(true);
		return 1.0e12;
	}
	if(sqr_r < _sqr_rcut) {
		return _patchy(p, q, compute_r, update_forces);
	}

	return 0.;
}

