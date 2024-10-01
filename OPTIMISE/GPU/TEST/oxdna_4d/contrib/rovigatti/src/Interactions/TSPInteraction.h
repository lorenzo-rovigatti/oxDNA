#ifndef TSP_INTERACTION_H
#define TSP_INTERACTION_H

#include <vector>
#include "Interactions/BaseInteraction.h"
#include "Particles/BaseParticle.h"
#include "Lists/Cells.h"

#define MAX_INSERTION_TRIES 1000000

class TSPParticle;

/**
 * @brief Handles interactions between TSPs.
 *
 * This interaction is selected with
 * interaction_type = TSP
 *
 * Input options: (type can take values 0, 1, 2 for A-A, A-B and B-B interactions, respectively)
 *
 * @verbatim
 TSP_rfene = <float> (FENE length constant for bonded interactions)
 TSP_sigma[type] = <float> (particle diameter associated to each interaction)
 TSP_epsilon[type] = <float> (energy scale associated to each interaction)
 TSP_attractive[type] = <float> (whether the interaction contains an attractive tail or not)
 TSP_n[type] = <int> (exponent for the generalised LJ potential for each interaction)
 [TSP_attractive_anchor = <bool> (set to true if you want the anchor monomer to be of type B instead of type A. Defaults to false)]
 [TSP_only_chains = <bool> (if true the system will be composed of chains only. The topology will be interpreted accordingly by ignoring the first value of each line (which, in the case of TSPs, is the number of arms). Defaults to false)]
 [TSP_only_intra = <bool> (if true monomers belonging to different stars will not interact. Defaults to false)]
 @endverbatim
 */
class TSPInteraction: public BaseInteraction {
protected:
	number _rfene, _sqr_rfene;
	number _rfene_anchor, _sqr_rfene_anchor;

	number _TSP_sqr_rep_rcut;
	number _TSP_lambda;
	int _TSP_n;

	number _TSP_yukawa_A;
	number _TSP_yukawa_xi;
	number _yukawa_E_cut;

	bool _attractive_anchor;
	bool _only_chains;
	bool _only_intra;
	bool _yukawa_repulsion;

	int *_N_arms;
	int *_N_monomer_per_arm;
	number *_alpha;
	int _N_stars;

	std::vector<TSPParticle *> _anchors;

	number _fene(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	number _nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	bool _insert_anchor(std::vector<BaseParticle *> &particles, BaseParticle *p, Cells *c);
	bool _does_overlap(std::vector<BaseParticle *> &particles, BaseParticle *p, Cells *c);

public:
	enum {
		BONDED = 0, NONBONDED = 1
	};
	TSPInteraction();
	virtual ~TSPInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	virtual int get_N_from_topology();
	virtual void read_topology(int *N_stars, std::vector<BaseParticle *> &particles);
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);

	virtual void generate_random_configuration(std::vector<BaseParticle *> &particles);

	void set_only_intra(bool value) {
		_only_intra = value;
	}
};

class TSPParticle: public BaseParticle {
protected:
	bool _is_anchor;
	int _arm;

public:
	TSPParticle();
	virtual ~TSPParticle();

	virtual bool is_rigid_body() { return false; }
	virtual bool is_bonded(BaseParticle *q);
	virtual bool is_anchor() { return _is_anchor; }
	virtual int arm() { return _arm; }

	virtual void add_bonded_neigh(TSPParticle *nn);
	virtual void flag_as_anchor() { _is_anchor = true; }
	virtual void set_arm(int na) { _arm = na; }

	std::set<TSPParticle *> bonded_neighs;
};

extern "C" TSPInteraction *make_TSPInteraction();

#endif /* TSP_INTERACTION_H */
