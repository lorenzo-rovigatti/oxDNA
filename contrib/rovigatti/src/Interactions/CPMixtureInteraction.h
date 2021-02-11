/*
 * CPMixtureInteraction.h
 *
 *  Created on: 10/feb/2013
 *      Author: lorenzo
 */

#ifndef CPMIXTUREINTERACTION_H_
#define CPMIXTUREINTERACTION_H_

#include "Interactions/BaseInteraction.h"

/**
 *
 */
class CPMixtureInteraction: public BaseInteraction {
protected:
	number _E_cut[3];
	int _CP_int_type[3];
	int _N_A, _N_B;
	int _n[3];

	number _sigma[3];
	number _sqr_sigma[3];
	number _epsilon[3];
	number _sqr_CP_rcut[3];
	number _CP_rcut[3];
	number _A[3], _B[3], _C[3], _D[3];

public:
	enum {
		CP = 0
	};

	enum {
		POWER_LAW = 0, GAUSSIAN = 1, LOUIS
	};

	CPMixtureInteraction();
	virtual ~CPMixtureInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);
	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
};

extern "C" CPMixtureInteraction *make_CPMixtureInteraction();

#endif /* CPMIXTUREINTERACTION_H_ */
