/*
 * FSInteraction.h
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#ifndef FSINTERACTION_H_
#define FSINTERACTION_H_

#include "BaseInteraction.h"

template <typename number>
struct FSBond {
	BaseParticle<number> *other;
	number r_p;
	int p_patch, q_patch;
	number energy;
	LR_vector<number> force;
	LR_vector<number> p_torque, q_torque;

	FSBond(BaseParticle<number> *o, number my_r_p, int pp, int qp, number e) : other(o), r_p(my_r_p), p_patch(pp), q_patch(qp), energy(e) {}
};

/**
 * @brief Manages the interaction between FS-like patchy particles
 *
 * This interaction is selected with
 * interaction_type = patchy
 *
 * @verbatim
FS_N = <int> (number of patches)
[FS_N_B = <int> (number of patches on species B)]
@endverbatim
 */
template <typename number>
class FSInteraction: public BaseInteraction<number, FSInteraction<number> > {
protected:
	/// Number of patches per particle
	int _N_patches;

	/// Number of patches per second-species particle
	int _N_patches_B;

	/// Number of particles per species
	int _N_A, _N_B, _N;

	/// true if we are to simulate a patchy binary mixture, false otherwise. If false then A-A interactions are enabled
	bool _one_component;

	/// Repulsive interaction energy at the cut-off
	number _rep_E_cut;
	number _rep_rcut;
	number _sqr_rep_rcut;

	number _patch_rcut;
	number _sqr_patch_rcut;
	number _sigma_ss;
	number _rcut_ss;
	number _lambda;
	number _A_part, _B_part;

	std::vector<std::vector<FSBond<number> > > _bonds;

	/**
	 * @brief FS interaction between two particles.
	 *
	 * @param p
	 * @param q
	 * @param r
	 * @param update_forces
	 * @return
	 */
	inline number _two_body(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	inline number _three_body(BaseParticle<number> *p, FSBond<number> &new_bond, bool update_forces);

public:
	enum {
		FS = 0
	};

	FSInteraction();
	virtual ~FSInteraction();

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

template<typename number>
number FSInteraction<number>::_two_body(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number sqr_r = r->norm();
	if(sqr_r > this->_sqr_rcut) return (number) 0.f;

	number energy = (number) 0.f;

	// centre-centre
	if(sqr_r < _sqr_rep_rcut) {
		number ir2 = 1. / sqr_r;
		number lj_part = ir2*ir2*ir2;
		energy = 4 * (SQR(lj_part) - lj_part) + _rep_E_cut;
		if(update_forces) {
			LR_vector<number> force = *r * (-24. * (lj_part - 2*SQR(lj_part)) / sqr_r);
			p->force -= force;
			q->force += force;
		}
	}

	if(p->type != q->type || _one_component) {
		LR_vector<number> tmptorquep(0, 0, 0);
		LR_vector<number> tmptorqueq(0, 0, 0);
		for(int pi = 0; pi < p->N_int_centers; pi++) {
			LR_vector<number> ppatch = p->int_centers[pi];
			for(int pj = 0; pj < q->N_int_centers; pj++) {
				LR_vector<number> qpatch = q->int_centers[pj];

				LR_vector<number> patch_dist = *r + qpatch - ppatch;
				number dist = patch_dist.norm();
				if(dist < _sqr_patch_rcut) {
					number r_p = sqrt(dist);
					number exp_part = exp(_sigma_ss / (r_p - _rcut_ss));
					number tmp_energy = _A_part * exp_part * (_B_part/SQR(dist) - 1.);

					energy += tmp_energy;

					FSBond<number> p_bond(q, r_p, pi, pj, tmp_energy);
					FSBond<number> q_bond(p, r_p, pj, pi, tmp_energy);

					if(update_forces) {
						number force_mod  = _A_part * exp_part * (4.*_B_part/(SQR(dist)*r_p)) + _sigma_ss * tmp_energy / SQR(r_p - _rcut_ss);
						LR_vector<number> tmp_force = patch_dist * (force_mod / r_p);

						LR_vector<number> p_torque = p->orientationT * ppatch.cross(tmp_force);
						LR_vector<number> q_torque = q->orientationT * qpatch.cross(tmp_force);

						p->force -= tmp_force;
						q->force += tmp_force;

						p->torque -= p_torque;
						q->torque += q_torque;

						p_bond.force = tmp_force;
						p_bond.p_torque = p_torque;
						p_bond.q_torque = q_torque;

						q_bond.force = -tmp_force;
						q_bond.p_torque = -q_torque;
						q_bond.q_torque = -p_torque;

					}

					energy += _three_body(p, p_bond, update_forces);
					energy += _three_body(q, q_bond, update_forces);

					_bonds[p->index].push_back(p_bond);
					_bonds[q->index].push_back(q_bond);
				}
			}
		}
	}

	return energy;
}

template<typename number>
number FSInteraction<number>::_three_body(BaseParticle<number> *p, FSBond<number> &new_bond, bool update_forces) {
	number energy = 0.;

	typename std::vector<FSBond<number> >::iterator it = _bonds[p->index].begin();
	for(; it != _bonds[p->index].end(); it++) {
		// three-body interactions happen only when the same patch is involved in
		// more than a bond
		if(it->p_patch == new_bond.p_patch) {
			number curr_energy = -new_bond.energy;
			if(new_bond.r_p < _sigma_ss) curr_energy = 1.;

			number other_energy = -it->energy;
			if(it->r_p < _sigma_ss) other_energy = 1.;

			energy += _lambda * curr_energy * other_energy;

			if(update_forces) {
				if(new_bond.r_p > _sigma_ss) {
					BaseParticle<number> *other = new_bond.other;

					number factor = -_lambda*other_energy;

					p->force -= factor*new_bond.force;
					other->force += factor*new_bond.force;

					p->torque -= factor*new_bond.p_torque;
					other->torque += factor*new_bond.q_torque;
				}

				if(it->r_p > _sigma_ss) {
					BaseParticle<number> *other = it->other;

					number factor = -_lambda*curr_energy;

					p->force -= factor*it->force;
					other->force += factor*it->force;

					p->torque -= factor*it->p_torque;
					other->torque += factor*it->q_torque;
				}
			}
		}
	}

	return energy;
}

extern "C" FSInteraction<float> *make_float() { return new FSInteraction<float>(); }
extern "C" FSInteraction<double> *make_double() { return new FSInteraction<double>(); }

#endif /* FSINTERACTION_H_ */
