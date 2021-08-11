/*
 * ValerioInteraction.h
 *
 *  Created on: 09/aug/2021
 *      Author: lorenzo
 */

#ifndef VALERIOINTERACTION_H_
#define VALERIOINTERACTION_H_

#include "Interactions/BaseInteraction.h"
#include "Particles/PatchyParticle.h"

class ValerioBond;

/**
 * @brief Manages the interaction between generic patchy particles.
 *
 * This interaction is selected with
 * interaction_type = ValerioInteraction
 */
class ValerioInteraction: public BaseInteraction {
protected:
	int _N = 0;

	number _E_cut;
	/// Patchy-related quantities
	number _sqr_patch_rcut;
	number _patch_alpha = 0.11;
	int _patch_pow = 10;
	number _patch_pow_alpha;

	/// 3-body semiflexibility
	number _3b_k = 10.;

	const int _rep_power = 200;

	number _patchy_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	number _spherical_patchy_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	std::set<ParticlePair> _bonds;

public:
	enum {
		PATCHY = 0,
		SPHERICAL = 1
	};

	ValerioInteraction();
	virtual ~ValerioInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	void begin_energy_computation() override;

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
};

extern "C" ValerioInteraction *make_ValerioInteraction();

#endif /* VALERIOINTERACTION_H_ */
