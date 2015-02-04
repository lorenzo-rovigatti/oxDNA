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
	number _centre_patch_dist;

	inline number _patchy(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
public:
	enum {
		HSC, PATCHY
	};

	PatchyHSC();
	virtual ~PatchyHSC();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
};

template<typename number>
number PatchyHSC<number>::_patchy(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number energy = 0.;
	for(int ppi = -1; ppi < 2; ppi += 2) {
		LR_vector<number> pp = p->pos + p->orientation.v3*(ppi*_centre_patch_dist);
		for(int pqi = -1; pqi < 2; pqi += 2) {
			LR_vector<number> pq = q->pos + q->orientation.v3*(pqi*_centre_patch_dist);

			if(pp.sqr_min_image_distance(pq, this->_box_side) < _sqr_patch_rcut) energy += -1.;
		}
	}

	return energy;
}

extern "C" PatchyHSC<float> *make_float() { return new PatchyHSC<float>(); }
extern "C" PatchyHSC<double> *make_double() { return new PatchyHSC<double>(); }

#endif /* PATCHYHSC_H_ */
