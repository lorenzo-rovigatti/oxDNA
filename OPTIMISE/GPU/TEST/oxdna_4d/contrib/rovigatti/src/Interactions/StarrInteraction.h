/*
 * StarrInteraction.h
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#ifndef STARRINTERACTION_H_
#define STARRINTERACTION_H_

#include "Interactions/BaseInteraction.h"

#include "Lists/Cells.h"

#define P_HUB 400

/**
 * @brief Manages the interaction between Starr tetramers.
 *
 * This interaction is selected with
 * interaction_type = StarrInteraction
 */
class StarrInteraction: public BaseInteraction {
protected:
	int _N_per_strand;
	int _N_strands, _N_tetramers, _N_dimers;
	int _N_dimer_spacers;

	bool _starr_model;
	int _mode;

	number _LJ_sigma[3];
	number _LJ_sqr_sigma[3];
	number _LJ_rcut[3];
	number _LJ_sqr_rcut[3];
	number _LJ_E_cut[3];
	number _der_LJ_E_cut[3];

	number _fene_K;
	number _fene_sqr_r0;

	number _lin_k;

	virtual number _fene(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	virtual number _three_body(BaseParticle *p, BaseParticle *n3, BaseParticle *n5, bool update_forces);

	virtual void _read_strand_topology(int *N_strands, std::vector<BaseParticle *> &particles);
	virtual void _read_tetramer_topology(int *N_strands, std::vector<BaseParticle *> &particles);
	virtual void _read_vitrimer_topology(int *N_strands, std::vector<BaseParticle *> &particles);

	virtual void _generate_strands(std::vector<BaseParticle *> &particles, Cells &c);
	virtual void _generate_tetramers(std::vector<BaseParticle *> &particles, Cells &c);
	virtual void _generate_vitrimers(std::vector<BaseParticle *> &particles, Cells &c);

public:
	enum {
		BONDED = 0, NONBONDED = 4
	};

	enum {
		STRANDS = 0, TETRAMERS = 1, VITRIMERS = 2
	};

	StarrInteraction();
	virtual ~StarrInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);

	virtual void generate_random_configuration(std::vector<BaseParticle *> &particles);
};

extern "C" StarrInteraction *make_StarrInteraction();

#endif /* STARRINTERACTION */
