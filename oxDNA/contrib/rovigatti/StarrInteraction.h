/*
 * StarrInteraction.h
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#ifndef STARRINTERACTION_H_
#define STARRINTERACTION_H_

#include "BaseInteraction.h"

#include "../Lists/Cells.h"

#define P_HUB 400

/**
 * @brief Manages the interaction between Starr tetramers.
 *
 * This interaction is selected with
 * interaction_type = StarrInteraction
 */
template <typename number>
class StarrInteraction: public BaseInteraction<number, StarrInteraction<number> > {
protected:
	int _N_per_strand;
	int _N_strands, _N_tetramers, _N_dimers;

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

	virtual number _fene(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _two_body(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _three_body(BaseParticle<number> *p, BaseParticle<number> *n3, BaseParticle<number> *n5, bool update_forces);

	virtual void _read_strand_topology(int N, int *N_strands, BaseParticle<number> **particles);
	virtual void _read_tetramer_topology(int N, int *N_strands, BaseParticle<number> **particles);
	virtual void _read_vitrimer_topology(int N, int *N_strands, BaseParticle<number> **particles);

	virtual void _generate_strands(BaseParticle<number> **particles, int N, number box_side, Cells<number> &c);
	virtual void _generate_tetramers(BaseParticle<number> **particles, int N, number box_side, Cells<number> &c);
	virtual void _generate_vitrimers(BaseParticle<number> **particles, int N, number box_side, Cells<number> &c);

public:
	enum {
		BONDED = 0,
		NONBONDED = 4
	};

	enum {
		STRANDS = 0,
		TETRAMERS = 1,
		VITRIMERS = 2
	};

	StarrInteraction();
	virtual ~StarrInteraction();

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

	virtual void generate_random_configuration(BaseParticle<number> **particles, int N, number box_side);
};

extern "C" StarrInteraction<float> *make_StarrInteraction_float();
extern "C" StarrInteraction<double> *make_StarrInteraction_double();

#endif /* STARRINTERACTION */
