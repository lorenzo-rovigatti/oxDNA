/*
 * StarrInteraction.h
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#ifndef STARRINTERACTION_H_
#define STARRINTERACTION_H_

#include "BaseInteraction.h"

/**
 * @brief Manages the interaction between FS-like patchy particles.
 *
 * This interaction is selected with
 * interaction_type = patchy
 *
 * @verbatim
FS_N = <int> (number of patches)
[FS_N_B = <int> (number of patches on species B)]
[FS_one_component = <bool> (if true then the system will not contain particles of type B and A-A bonds will be allowed. Defaults to false.)]
@endverbatim
 */
template <typename number>
class StarrInteraction: public BaseInteraction<number, StarrInteraction<number> > {
protected:
	int _N_per_tetramer;
	int _N_per_strand;
	int _N_strands_per_tetramer;
	int _N_tetramers;

	number _LJ_sigma[3];
	number _LJ_sqr_sigma[3];
	number _LJ_sqr_rcut[3];
	number _LJ_E_cut[3];
	number _der_LJ_E_cut[3];

	number _fene_K;
	number _fene_sqr_r0;

	number _lin_k;

	virtual number _fene(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _two_body(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _three_body(BaseParticle<number> *p, bool update_forces);

public:
	enum {
		STARR = 0
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

//	virtual void generate_random_configuration(BaseParticle<number> **particles, int N, number box_side);
};

extern "C" StarrInteraction<float> *make_StarrInteraction_float();
extern "C" StarrInteraction<double> *make_StarrInteraction_double();

#endif /* STARRINTERACTION */
