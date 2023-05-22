/*
 * PatchyParticle.h
 *
 *  Created on: 15/mar/2013
 *      Author: lorenzo
 */

#ifndef PATCHYPARTICLE_H_
#define PATCHYPARTICLE_H_

#include "BaseParticle.h"

struct PatchyBond {
	BaseParticle *other;
	number r_p;
	int p_patch, q_patch;
	number energy;
	LR_vector force;
	LR_vector p_torque, q_torque;

	PatchyBond(BaseParticle *o, number my_r_p, int pp, int qp, number e) :
					other(o),
					r_p(my_r_p),
					p_patch(pp),
					q_patch(qp),
					energy(e) {
	}
};

/**
 * @brief Incapsulates a patchy particle.
 */

class PatchyParticle : public BaseParticle {
protected:
	number _sigma;
	std::vector<LR_vector> _base_patches;
	int _N_patches;

public:
	PatchyParticle(std::vector<LR_vector> base_patches, int nt, number sigma);
	PatchyParticle(int N_patches, int nt, number sigma);
	virtual ~PatchyParticle();

	void set_positions();

	virtual bool is_rigid_body() {
		return true;
	}

	uint N_int_centers() override {
		return _N_patches;
	}

	std::vector<LR_vector> base_patches() {
		return _base_patches;
	}

	std::vector<PatchyBond> bonds;
};

#endif /* PATCHYPARTICLE_H_ */
