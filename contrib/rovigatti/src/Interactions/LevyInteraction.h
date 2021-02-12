/*
 * LevyInteraction.h
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#ifndef LEVYINTERACTION_H_
#define LEVYINTERACTION_H_

#include "Interactions/BaseInteraction.h"

#include "Lists/Cells.h"
#include "Particles/CustomParticle.h"

/**
 * @brief Manages the interaction between Starr tetramers.
 *
 * This interaction is selected with
 * interaction_type = LevyInteraction
 */
class LevyInteraction: public BaseInteraction {
protected:
	int _N_tetramers, _N_dimers, _N_monomers;
	int _N_per_tetramer, _N_per_dimer;
	int _patchy_power;
	bool _rigid_model;

	number _sigma[3];
	number _sqr_sigma[3];
	number _sqr_tot_rcut[3];
	number _E_cut[3];
	number _epsilon, _monomer_epsilon;
	number _patch_E_cut, _patch_monomer_E_cut;
	number _sqr_patch_rcut;
	number _patch_alpha;
	number _patch_pow_alpha;

	number _fene_K;
	number _fene_sqr_r0;

	number _lin_k, _terminal_lin_k;

	virtual number _fene(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _three_body(BaseParticle *p, BaseParticle *n3, BaseParticle *n5, bool update_forces);

public:
	enum {
		BONDED = 0, NONBONDED = 4
	};

	enum {
		TETRA_CENTRE = 0, TETRA_PATCHY = 1, DIMER_CENTRE = 2, DIMER_PATCHY = 4, MONOMER = 8
	};

	LevyInteraction();
	virtual ~LevyInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);

//	virtual void generate_random_configuration(std::vector<BaseParticle *> &particles, number box_side);
};

class PatchySite: public CustomParticle {
protected:
	number _sigma, _radius;

public:
	PatchySite() {
		this->int_centers.resize(1);
		_sigma = 1.;
		_radius = _sigma * 0.5;
	}

	virtual ~PatchySite() {

	}

	virtual bool is_rigid_body() {
		return true;
	}
	virtual void set_sigma(number s) {
		_sigma = s;
		_radius = s * 0.5;
	}

	void set_positions() {
		// this is equivalent to orientation*(1, 0, 0)
		this->int_centers[0] = this->orientationT.v1 * _radius;
	}
};

extern "C" LevyInteraction *make_LevyInteraction();

#endif /* LEVYINTERACTION */
