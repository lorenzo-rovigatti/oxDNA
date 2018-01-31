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
template <typename number>
class LevyInteraction: public BaseInteraction<number, LevyInteraction<number> > {
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

	virtual number _fene(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _two_body(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _three_body(BaseParticle<number> *p, BaseParticle<number> *n3, BaseParticle<number> *n5, bool update_forces);

public:
	enum {
		BONDED = 0,
		NONBONDED = 4
	};

	enum {
		TETRA_CENTRE = 0,
		TETRA_PATCHY = 1,
		DIMER_CENTRE = 2,
		DIMER_PATCHY = 4,
		MONOMER = 8
	};

	LevyInteraction();
	virtual ~LevyInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(BaseParticle<number> **particles, int N);

	virtual number pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_term(int name, BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, r, update_forces);
	}

	virtual void read_topology(int N, int *N_strands, BaseParticle<number> **particles);
	virtual void check_input_sanity(BaseParticle<number> **particles, int N);

//	virtual void generate_random_configuration(BaseParticle<number> **particles, int N, number box_side);
};

template<typename number>
class PatchySite: public CustomParticle<number> {
protected:
	number _sigma, _radius;

public:
	PatchySite() {
		this->N_int_centers = 1;
		this->int_centers = new LR_vector<number>[1];
		_sigma = 1.;
		_radius = _sigma*0.5;
	}

	virtual ~PatchySite() {

	}

	virtual bool is_rigid_body() { return true; }
	virtual void set_sigma(number s) {
		_sigma = s;
		_radius = s*0.5;
	}

	void set_positions() {
		// this is equivalent to orientation*(1, 0, 0)
		this->int_centers[0] = this->orientationT.v1*_radius;
	}
};

extern "C" LevyInteraction<float> *make_LevyInteraction_float();
extern "C" LevyInteraction<double> *make_LevyInteraction_double();

#endif /* LEVYINTERACTION */
