#ifndef DETAILEDPOLYMERSWAP_INTERACTION_H
#define DETAILEDPOLYMERSWAP_INTERACTION_H

#include "Interactions/BaseInteraction.h"

#include <vector>
#include <array>

struct PSBond {
	BaseParticle *other;
	number energy;
	number epsilon;
	LR_vector force;
	LR_vector r;

	PSBond(BaseParticle *o, number e, number eps, LR_vector nr) :
		other(o),
		energy(e),
		epsilon(eps),
		r(nr) {
			
	}
};

struct PSBondCompare {
	bool operator()(const PSBond &lhs, const PSBond &rhs) {
		if(lhs.other->index == rhs.other->index) {
			return false;
		}
		else {
			return true;
		}
	}
};

/**
 * @brief Handles interactions in microgel systems.
 *
 * This interaction is selected with
 * interaction_type = DetailedPolymerSwapInteraction
 *
 * @verbatim

 @endverbatim
 */
class DetailedPolymerSwapInteraction: public BaseInteraction {
protected:
	number _Kfene = 15.;
	number _rfene = 1.5;
	number _sqr_rfene;
	number _WCA_sigma = 1.0;
	number _PS_sqr_rep_rcut;

	std::vector<LR_vector> _chain_coms;

	number _PS_alpha = 0.;
	number _PS_beta = 0.;
	number _PS_gamma = 0.;
	int _PS_n = 6;

	int _N_attractive_types = 0;
	int _interaction_matrix_size = 0;
	std::string _interaction_matrix_file;

	/// three-body potential stuff
	number _3b_rcut = -1.;
	number _sqr_3b_rcut = -1.;
	number _3b_sigma = 1.05;
	number _3b_range = 1.6;
	number _3b_lambda = 1.0;
	std::vector<number> _3b_epsilon;
	number _3b_prefactor = 0.1;
	number _3b_A_part = 0.;
	number _3b_B_part = 0.;

	/// three-body flexibility stuff
	bool _enable_semiflexibility = false;
	number _semiflexibility_k;

	std::map<int, std::set<PSBond, PSBondCompare> > _bonds;

	int _N_chains = -1;

	StressTensor _inter_chain_stress_tensor;

	void _parse_interaction_matrix();

	void _update_inter_chain_stress_tensor(int chain, int ref_chain, LR_vector group_force);

	bool _sticky_interaction(int p_btype, int q_btype);
	number _fene(BaseParticle *p, BaseParticle *q, bool update_forces);
	number _WCA(BaseParticle *p, BaseParticle *q, bool update_forces);
	number _sticky(BaseParticle *p, BaseParticle *q, bool update_forces);
	number _patchy_three_body(BaseParticle *p, PSBond &new_bond, bool update_forces);
	number _semiflexibility_three_body(BaseParticle *middle, BaseParticle *n1, BaseParticle *n2, bool update_forces);

public:
	enum {
		BONDED = 0, NONBONDED = 1, STICKY = 2
	};
	enum {
		MONOMER = 0, STICKY_ANY = 1
	};

	bool no_three_body = false;

	DetailedPolymerSwapInteraction();
	virtual ~DetailedPolymerSwapInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	number P_inter_chain();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	void begin_energy_computation() override;
	void begin_energy_and_force_computation() override;

	bool has_custom_stress_tensor() const override;

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	virtual number pair_nonbonded_WCA(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_nonbonded_sticky(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	virtual void read_topology(int *N_stars, std::vector<BaseParticle *> &particles);
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
};

extern "C" DetailedPolymerSwapInteraction *make_DetailedPolymerSwapInteraction();

#endif /* DETAILEDPOLYMERSWAP_INTERACTION_H */
