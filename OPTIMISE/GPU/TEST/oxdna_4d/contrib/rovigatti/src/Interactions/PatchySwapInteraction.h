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

/**
 * @brief Manages the interaction between generic patchy particles.
 *
 * This interaction is selected with
 * interaction_type = PatchySwapInteraction
 */
class PatchySwapInteraction: public BaseInteraction {
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

	void _parse_interaction_matrix();
	std::vector<LR_vector> _parse_base_patches(std::string filename, int N_patches);

	number _patchy_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	number _spherical_patchy_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	number _three_body(BaseParticle *p, PatchyBond &new_bond, bool update_forces);

	inline std::vector<PatchyBond> &_particle_bonds(BaseParticle *p) {
		return static_cast<PatchyParticle *>(p)->bonds;
	}

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

	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
};

extern "C" PatchySwapInteraction *make_PatchySwapInteraction();

#endif /* PATCHYSWAPINTERACTION_H_ */
