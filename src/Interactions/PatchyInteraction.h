/*
 * PatchyInteraction.h
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#ifndef PATCHYINTERACTION_H_
#define PATCHYINTERACTION_H_

#define PATCHY_POWER 200

#include "BaseInteraction.h"
#include "../Particles/PatchyParticle.h"

/**
 * @brief Manages the interaction between simple patchy particles (as described in http://jcp.aip.org/resource/1/jcpsa6/v131/i1/p014504_s1)
 *
 * This interaction is selected with
 * interaction_type = patchy
 *
 * @verbatim
PATCHY_N = <int> (number of patches)
[PATCHY_N_B = <int> (number of patches on species B)]
[PATCHY_alpha = <float> (width of patches, defaults to 0.12)]
[PATCHY_epsilon_AA = <float> (depth of the well of the patch-patch interaction between particles of species A)]
[PATCHY_epsilon_BB = <float> (depth of the well of the patch-patch interaction between particles of species B)]
[PATCHY_epsilon_AB = <float> (depth of the well of the patch-patch interaction between particles of different species)]
[PATCHY_sigma_AA = <float> (diameter controlling the repulsive interaction between particles of species A)]
[PATCHY_sigma_BB = <float> (diameter controlling the repulsive interaction between particles of species B)]
[PATCHY_sigma_AB = <float> (diameter controlling the repulsive interaction between particles of different species)]
@endverbatim
 */
class PatchyInteraction: public BaseInteraction {
protected:
	/// Number of patches per particle
	int _N_patches;

	/// Number of patches per second-species particle
	int _N_patches_B;

	/// Number of particles per species
	int _N_A, _N_B;

	/// True if we are to simulate a patchy binary mixture, false otherwise
	bool _is_binary;

	/// Particles' diameters
	number _sigma[3];

	/// Squared diameters
	number _sqr_sigma[3];

	/// Repulsive interaction cut-off squared
	number _sqr_tot_rcut[3];

	/// Depth of the patch-patch well
	number _epsilon[3];

	/// Repulsive interaction energy at the cut-off
	number _E_cut[3];

	/// Patch-patch interaction energy at the cut-off
	number _patch_E_cut[3];

	/// Patch-patch interaction cut-off squared
	number _sqr_patch_rcut;

	/// Width of the patch, defaults to 0.12
	number _patch_alpha;

	/// _patch_alpha^10
	number _patch_pow_alpha;

	number _bond_energy_threshold = -0.1;

	/**
	 * @brief Patchy interaction between two particles.
	 *
	 * @param p
	 * @param q
	 * @param r
	 * @param update_forces
	 * @return
	 */
	number _patchy_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	inline std::vector<PatchyBond> &_particle_bonds(BaseParticle *p) {
		return static_cast<PatchyParticle *>(p)->bonds;
	}

public:
	enum {
		PATCHY = 4
	};

	PatchyInteraction();
	virtual ~PatchyInteraction();

	void get_settings(input_file &inp) override;
	void init() override;

	number get_alpha() { return _patch_alpha; }

	void allocate_particles(std::vector<BaseParticle *> &particles) override;

	void begin_energy_computation() override;

	number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false) override;
	number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false) override;
	number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces=false) override;

	void read_topology(int *N_strands, std::vector<BaseParticle *> &particles) override;
	void check_input_sanity(std::vector<BaseParticle *> &particles) override;
};

#endif /* PATCHYINTERACTION_H_ */
