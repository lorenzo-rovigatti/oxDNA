/*
 * PatchyInteraction.h
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#ifndef PATCHYINTERACTION_H_
#define PATCHYINTERACTION_H_

#define PATCHY_POWER 200
#define PATCHY_CUTOFF 0.18f

#include "BaseInteraction.h"

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
@endverbatim
 */
template <typename number>
class PatchyInteraction: public BaseInteraction<number, PatchyInteraction<number> > {
protected:
	/// Number of patches per particle
	int _N_patches;

	/// Number of patches per second-species particle
	int _N_patches_B;

	/// Number of particles per species
	int _N_A, _N_B;

	/// true if we are to simulate a patchy binary mixture, false otherwise
	bool _is_binary;

	/// Repulsive interaction energy at the cut-off
	number _E_cut;

	/// Patch-patch interaction energy at the cut-off
	number _patch_E_cut;

	/// Width of the patch, defaults to 0.12
	number _patch_alpha;

	/// _patch_alpha^10
	number _patch_pow_alpha;

	/**
	 * @brief Patchy interaction between two particles.
	 *
	 * @param p
	 * @param q
	 * @param r
	 * @param update_forces
	 * @return
	 */
	inline number _patchy_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);

public:
	enum {
		PATCHY = 0
	};

	PatchyInteraction();
	virtual ~PatchyInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	number get_alpha() { return _patch_alpha; }

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

template<typename number>
number PatchyInteraction<number>::_patchy_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number rnorm = r->norm();
	if(rnorm > this->_sqr_rcut) return (number) 0.f;

	number energy = (number) 0.f;

	number part = 1 / powf(rnorm, PATCHY_POWER * 0.5f);
	energy = part - _E_cut;

	if(update_forces) {
		LR_vector<number> force = *r * (PATCHY_POWER * part / rnorm);
		p->force -= force;
		q->force += force;
	}

	int c = 0;
	LR_vector<number> tmptorquep(0, 0, 0);
	LR_vector<number> tmptorqueq(0, 0, 0);
	for(int pi = 0; pi < p->N_int_centers; pi++) {
		LR_vector<number> ppatch = p->int_centers[pi];
		for(int pj = 0; pj < q->N_int_centers; pj++) {
			LR_vector<number> qpatch = q->int_centers[pj];

			LR_vector<number> patch_dist = *r + qpatch - ppatch;
			number dist = patch_dist.norm();
			if(dist < SQR(PATCHY_CUTOFF)) {
				c++;
				number r8b10 = dist*dist*dist*dist / _patch_pow_alpha;
				number exp_part = -1.001f * exp(-(number)0.5f * r8b10 * dist);

				energy += exp_part - _patch_E_cut;

				if(update_forces) {
					LR_vector<number> tmp_force = patch_dist * (5 * exp_part * r8b10);

					p->torque -= p->orientationT * ppatch.cross(tmp_force);
					q->torque += q->orientationT * qpatch.cross(tmp_force);

					p->force -= tmp_force;
					q->force += tmp_force;

					// we are in the single bond per patch condition and hence we can safely return here
					return energy;
				}
			}
		}
	}

	return energy;
}

#endif /* PATCHYINTERACTION_H_ */
