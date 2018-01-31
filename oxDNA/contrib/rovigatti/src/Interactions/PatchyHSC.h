/*
 * PatchyHSC.h
 *
 *  Created on: 11/ago/2014
 *      Author: lorenzo
 */

#ifndef PATCHYHSC_H_
#define PATCHYHSC_H_

#include "HardSpheroCylinderInteraction.h"

template <typename number>
class PatchyHSC: public HardSpheroCylinderInteraction<number> {
protected:
	number _protrusion;
	number _patch_r;
	number _sqr_patch_rcut;
	number _sqr_patch_shoulder_rcut;
	number _centre_patch_dist;
	number _spherocylinder_length, _sqr_spherocylinder_length;

	inline number _patchy_energy(LR_vector<number> r);
	inline number _patchy(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
public:
	enum {
		HSC, PATCHY
	};

	PatchyHSC();
	virtual ~PatchyHSC();

	virtual void get_settings(input_file &inp);
	virtual void init();
	virtual void allocate_particles(BaseParticle<number> **particles, int N);

	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
};

template<typename number>
number PatchyHSC<number>::_patchy_energy(LR_vector<number> r) {
	number dist = r.norm();
	if(dist < _sqr_patch_shoulder_rcut) return -0.5;
	else if(dist < _sqr_patch_rcut) return -1.0;
	return 0.0;
}

template<typename number>
number PatchyHSC<number>::_patchy(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number energy = 0.;
	// here we try pp + pq and -(pp + pq)
	LR_vector<number> ppq = p->int_centers[0] + q->int_centers[0];
	energy += _patchy_energy(*r + ppq);
	energy += _patchy_energy(*r - ppq);
	if(energy < 0.) return energy;

	// and here pp - pq and -(pp - pq)
	LR_vector<number> pmq = p->int_centers[0] - q->int_centers[0];
	energy += _patchy_energy(*r + pmq);
	energy += _patchy_energy(*r - pmq);

	return energy;
}

template<typename number>
class PatchySpherocylinder: public BaseParticle<number> {
private:
	number _centre_patch_dist;

public:
	PatchySpherocylinder(number cp_dist) : BaseParticle<number>(), _centre_patch_dist(cp_dist) {
		this->N_int_centers = 1;
		this->int_centers = new LR_vector<number>[1];
	};

	virtual ~PatchySpherocylinder() {};

	virtual bool is_rigid_body() { return true; }
	virtual void set_positions() {
		this->int_centers[0] = this->orientation.v3*_centre_patch_dist;
	}
};

extern "C" PatchyHSC<float> *make_float() { return new PatchyHSC<float>(); }
extern "C" PatchyHSC<double> *make_double() { return new PatchyHSC<double>(); }

#endif /* PATCHYHSC_H_ */
