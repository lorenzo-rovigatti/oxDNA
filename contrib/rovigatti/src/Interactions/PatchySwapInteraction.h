/*
 * PatchySwapInteraction.h
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#ifndef PATCHYSWAPINTERACTION_H_
#define PATCHYSWAPINTERACTION_H_

#include "Interactions/BaseInteraction.h"
#include "Particles/PatchyParticle.h"

struct FSBond {
	BaseParticle *other;
	number r_p;
	int p_patch, q_patch;
	number energy;
	LR_vector force;
	LR_vector p_torque, q_torque;

	FSBond(BaseParticle *o, number my_r_p, int pp, int qp, number e) :
					other(o),
					r_p(my_r_p),
					p_patch(pp),
					q_patch(qp),
					energy(e) {
	}
};

struct FSBondCompare {
	bool operator()(const FSBond &lhs, const FSBond &rhs) {
		if(lhs.other->index == rhs.other->index && lhs.p_patch == rhs.p_patch && lhs.q_patch == rhs.q_patch) return false;
		else return true;
	}
};

/**
 * @brief Manages the interaction between FS-like patchy particles.
 *
 * This interaction is selected with
 * interaction_type = PatchySwapInteraction
 *
 * @verbatim
 FS_N = <int> (number of patches)
 [FS_N_B = <int> (number of patches on species B)]
 [FS_one_component = <bool> (if true then the system will not contain particles of type B and A-A bonds will be allowed. Defaults to false.)]
 @endverbatim
 */
class PatchySwapInteraction: public BaseInteraction<PatchySwapInteraction> {
protected:
	/// Number of particles of each species
	std::vector<int> _N_per_species;
	int _N_species = 0;
	int _N = 0;

	/// Number of patches per particle
	std::vector<int> _N_patches;
	std::vector<std::vector<LR_vector>> _base_patches;

	std::string _interaction_matrix_file;

	/// Repulsive interaction energy at the cut-off
	number _rep_rcut = 0.;
	number _sqr_rep_rcut = -1.;

	/// Patchy-related quantities
	std::vector<number> _patchy_eps;
	number _patch_rcut = -1.;
	number _sqr_patch_rcut = -1.;
	number _sigma_ss = 0.4;
	number _rcut_ss = -1.;
	number _lambda = 1.0;
	number _A_part = 0.;
	number _B_part = 0.;

	/// Optional spherical attraction
	number _spherical_attraction_strength = 0.;
	number _spherical_rcut = 2.5;
	number _sqr_spherical_rcut = 6.25;
	number _spherical_E_cut = 0.;

	std::vector<std::vector<FSBond>> _bonds;

	void _parse_interaction_matrix();
	std::vector<LR_vector> _parse_base_patches(std::string filename, int N_patches);

	number _patchy_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	number _spherical_patchy_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	number _three_body(BaseParticle *p, FSBond &new_bond, bool update_forces);

public:
	enum {
		PATCHY = 0,
		SPHERICAL = 1
	};

	bool no_three_body = false;

	PatchySwapInteraction();
	virtual ~PatchySwapInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	void begin_energy_computation() override;

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_term(int name, BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, compute_r, update_forces);
	}

	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
};

extern "C" PatchySwapInteraction *make_PatchySwapInteraction();

#endif /* PATCHYSWAPINTERACTION_H_ */
