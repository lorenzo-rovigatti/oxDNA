/*
 * mWInteraction.h
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#ifndef MWINTERACTION_H_
#define MWINTERACTION_H_

#include "Interactions/BaseInteraction.h"

struct mWBond {
	BaseParticle *other;
	LR_vector r;
	number mod_r;

	mWBond(BaseParticle *o, LR_vector my_r, number my_mod_r) :
					other(o),
					r(my_r),
					mod_r(my_mod_r) {
	}
};

/**
 * @brief Manages the interaction between mW-like patchy particles.
 *
 * This interaction is selected with
 * interaction_type = patchy
 *
 * @verbatim
 mW_N = <int> (number of patches)
 [mW_N_B = <int> (number of patches on species B)]
 [mW_one_component = <bool> (if true then the system will not contain particles of type B and A-A bonds will be allowed. Defaults to false.)]
 @endverbatim
 */
class mWInteraction: public BaseInteraction {
protected:
	int _N;

	number _A, _B;
	number _a;
	number _gamma;
	number _lambda;
	number _theta0, _cos_theta0;

	std::vector<std::vector<mWBond> > _bonds;

	/**
	 * @brief mW interaction between two particles.
	 *
	 * @param p
	 * @param q
	 * @param r
	 * @param update_forces
	 * @return
	 */
	inline number _two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	inline number _three_body(BaseParticle *p, mWBond &new_bond, bool update_forces);

public:
	enum {
		mW = 0
	};

	mWInteraction();
	virtual ~mWInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	void begin_energy_computation() override;

	virtual bool has_custom_stress_tensor() const {
		return true;
	}

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
};

extern "C" mWInteraction *make_mWInteraction();

#endif /* MWINTERACTION_H_ */
