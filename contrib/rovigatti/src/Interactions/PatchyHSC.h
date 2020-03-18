/*
 * PatchyHSC.h
 *
 *  Created on: 11/ago/2014
 *      Author: lorenzo
 */

#ifndef PATCHYHSC_H_
#define PATCHYHSC_H_

#include "HardSpheroCylinderInteraction.h"

class PatchyHSC: public HardSpheroCylinderInteraction {
protected:
	number _protrusion;
	number _patch_r;
	number _sqr_patch_rcut;
	number _sqr_patch_shoulder_rcut;
	number _centre_patch_dist;
	number _spherocylinder_length, _sqr_spherocylinder_length;

	inline number _patchy_energy(LR_vector r);
	inline number _patchy(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
public:
	enum {
		HSC, PATCHY
	};

	PatchyHSC();
	virtual ~PatchyHSC();

	virtual void get_settings(input_file &inp);
	virtual void init();
	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
};

number PatchyHSC::_patchy_energy(LR_vector r) {
	number dist = r.norm();
	if(dist < _sqr_patch_shoulder_rcut) return -0.5;
	else if(dist < _sqr_patch_rcut) return -1.0;
	return 0.0;
}

number PatchyHSC::_patchy(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = 0.;
	// here we try pp + pq and -(pp + pq)
	LR_vector ppq = p->int_centers[0] + q->int_centers[0];
	energy += _patchy_energy(_computed_r + ppq);
	energy += _patchy_energy(_computed_r - ppq);
	if(energy < 0.) return energy;

	// and here pp - pq and -(pp - pq)
	LR_vector pmq = p->int_centers[0] - q->int_centers[0];
	energy += _patchy_energy(_computed_r + pmq);
	energy += _patchy_energy(_computed_r - pmq);

	return energy;
}

class PatchySpherocylinder: public BaseParticle {
private:
	number _centre_patch_dist;

public:
	PatchySpherocylinder(number cp_dist) :
					BaseParticle(),
					_centre_patch_dist(cp_dist) {
		this->N_int_centers = 1;
		this->int_centers = new LR_vector[1];
	}
	;

	virtual ~PatchySpherocylinder() {
	}
	;

	virtual bool is_rigid_body() {
		return true;
	}
	virtual void set_positions() {
		this->int_centers[0] = this->orientation.v3 * _centre_patch_dist;
	}
};

extern "C" PatchyHSC *make_PatchyHSC() {
	return new PatchyHSC();
}

#endif /* PATCHYHSC_H_ */
