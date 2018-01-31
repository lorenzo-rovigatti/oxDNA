/*
 * FSInteraction.h
 *
 *  Created on: 14/mar/2013
 *      Author: lorenzo
 */

#ifndef FSINTERACTION_H_
#define FSINTERACTION_H_

#include "Interactions/BaseInteraction.h"

template <typename number>
struct FSBond {
	BaseParticle<number> *other;
	LR_vector<number> r;
	number r_p;
	int p_patch, q_patch;
	number energy;
	LR_vector<number> force;
	LR_vector<number> p_torque, q_torque;

	FSBond(BaseParticle<number> *o, LR_vector<number> my_r, number my_r_p, int pp, int qp, number e) : other(o), r(my_r), r_p(my_r_p), p_patch(pp), q_patch(qp), energy(e) {}
};

template <typename number>
struct FSBondCompare {
	bool operator() (const FSBond<number> &lhs, const FSBond<number> &rhs) {
		if(lhs.other->index == rhs.other->index && lhs.p_patch == rhs.p_patch && lhs.q_patch == rhs.q_patch) return false;
		else return true;
	}
};

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
class FSInteraction: public BaseInteraction<number, FSInteraction<number> > {
protected:
	/// Number of patches per particle
	int _N_patches;

	/// Number of patches per second-species particle
	int _N_patches_B;

	/// Number of particles per species
	int _N_A, _N_B, _N;

	/// Number of defective first-species particles
	int _N_def_A;

	/// true if we are to simulate a patchy binary mixture, false otherwise. If false then A-A interactions are enabled
	bool _one_component;

	/// true if particles of type B can interact between themselves
	bool _B_attraction;

	/// true if the two types of particle bear the same patches
	bool _same_patches;

	bool _needs_reset;

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

	std::vector<std::set<FSBond<number>, FSBondCompare<number> > > _bonds;

	std::vector<std::vector<number> > _stress_tensor;

	/**
	 * @brief FS interaction between two particles.
	 *
	 * @param p
	 * @param q
	 * @param r
	 * @param update_forces
	 * @return
	 */
	number _two_body(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);

	number _three_body(BaseParticle<number> *p, FSBond<number> &new_bond, bool update_forces);

	void _update_stress_tensor(LR_vector<number> r, LR_vector<number> f);

	bool _attraction_allowed(int p_type, int q_type);

public:
	enum {
		FS = 0
	};

	bool no_three_body;

	std::vector<std::vector<number> > get_stress_tensor() { return _stress_tensor; }

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
};

template<typename number>
void FSInteraction<number>::_update_stress_tensor(LR_vector<number> r, LR_vector<number> f) {
	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++) {
			number ri = r[i];
			number fj = f[j];
			_stress_tensor[i][j] += ri*fj;
		}
	}
}

template<typename number>
bool FSInteraction<number>::_attraction_allowed(int p_type, int q_type) {
	if(_same_patches) return true;
	if(_one_component) return true;
	if(p_type != q_type) return true;
	if(_B_attraction && p_type == P_B && q_type == P_B) return true;
	return false;
}

extern "C" FSInteraction<float> *make_FSInteraction_float();
extern "C" FSInteraction<double> *make_FSInteraction_double();

#endif /* FSINTERACTION_H_ */
