/*
 * FSInteraction.h
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#ifndef FSINTERACTION_H_
#define FSINTERACTION_H_

#include "Interactions/BaseInteraction.h"
#include "Particles/PatchyParticle.h"

/**
 * @brief Manages the interaction between FS-like patchy particles.
 *
 * This interaction is selected with
 * interaction_type = FSInteraction
 *
 * @verbatim
 FS_N = <int> (number of patches)
 [FS_N_B = <int> (number of patches on species B)]
 [FS_one_component = <bool> (if true then the system will not contain particles of type B and A-A bonds will be allowed. Defaults to false.)]
 @endverbatim
 */
class FSInteraction: public BaseInteraction {
protected:
	/// Number of patches per particle
	int _N_patches = -1;

	/// Number of patches per second-species particle
	int _N_patches_B = -1;

	/// Number of particles per species
	int _N_A = 0;
	int _N_B = 0;
	int _N = 0;

	/// Number of defective first-species particles
	int _N_def_A = 0;

	/// true if we are to simulate a patchy binary mixture, false otherwise. If false then A-A interactions are enabled
	bool _one_component = false;

	/// true if particles of type B can interact between themselves
	bool _B_attraction = false;

	/// true if the two types of particle bear the same patches
	bool _same_patches = false;

	/// Repulsive interaction energy at the cut-off
	number _rep_rcut = 0.;
	number _sqr_rep_rcut = -1.;

	/// Patchy-related quantities
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

	// used if we also have polymers in the simulation
	bool _with_polymers = false;
	int _N_in_polymers = 0;
	std::string _bond_filename;
	number _polymer_length_scale = 1.0;
	number _polymer_length_scale_sqr = 1.0;
	number _polymer_rfene = 1.5;
	number _polymer_rfene_sqr = 2.25;
	number _polymer_alpha = 0;
	number _polymer_beta = 0;
	number _polymer_gamma = 0;
	/// this quantity rescales the energy of the WCA and FENE parts. It is set to the thermal energy by default
	number _polymer_energy_scale = 1.0;

	/**
	 * @brief FS interaction between two particles.
	 *
	 * @param p
	 * @param q
	 * @param r
	 * @param update_forces
	 * @return
	 */
	number _patchy_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	number _spherical_patchy_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	number _polymer_fene(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	number _polymer_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	number _three_body(BaseParticle *p, PatchyBond &new_bond, bool update_forces);

	void _parse_bond_file(std::vector<BaseParticle *> &particles);

	bool _attraction_allowed(int p_type, int q_type);

	inline bool _is_patchy_patchy(int p_type, int q_type);
	inline bool _is_patchy_polymer(int p_type, int q_type);
	inline bool _is_polymer_polymer(int p_type, int q_type);

	inline std::vector<PatchyBond> &_particle_bonds(BaseParticle *p) {
		if(p->type != POLYMER) {
			return static_cast<PatchyParticle *>(p)->bonds;
		}
		else {
			throw oxDNAException("Cannot get a list of bonds from a particle that is not a PatchyParticle");
		}
	}

public:
	enum {
		PATCHY = 0,
		SPHERICAL = 1
	};

	enum {
		PATCHY_A = 0, PATCHY_B = 1, POLYMER = 2
	};

	bool no_three_body = false;

	FSInteraction();
	virtual ~FSInteraction();

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

extern "C" FSInteraction *make_FSInteraction();

#endif /* FSINTERACTION_H_ */
