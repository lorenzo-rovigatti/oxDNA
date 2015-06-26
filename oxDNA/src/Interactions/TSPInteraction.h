
#ifndef TSP_INTERACTION_H
#define TSP_INTERACTION_H

#include <vector>
#include "BaseInteraction.h"
#include "../Particles/TSPParticle.h"
#include "../Lists/Cells.h"

#define MAX_INSERTION_TRIES 1000000

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
template <typename number>
class TSPInteraction : public BaseInteraction<number, TSPInteraction<number> > {
protected:
	number _rfene, _sqr_rfene;
	number _rfene_anchor, _sqr_rfene_anchor;

	number _TSP_sqr_rep_rcut;
	number _TSP_lambda;
	int _TSP_n;

	bool _attractive_anchor;
	bool _only_chains;
	bool _only_intra;

	int *_N_arms;
	int *_N_monomer_per_arm;
	number *_alpha;
	int _N_stars;

	std::vector<TSPParticle<number> *> _anchors;

	number _fene(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	number _nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);

	bool _insert_anchor(BaseParticle<number> **particles, BaseParticle<number> *p, Cells<number> *c);
	bool _does_overlap(BaseParticle<number> **particles, BaseParticle<number> *p, Cells<number> *c);

public:
	enum {
		BONDED = 0,
		NONBONDED = 1
	};
	TSPInteraction();
	virtual ~TSPInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(BaseParticle<number> **particles, int N);

	virtual number pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_term(int name, BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, r, update_forces);
	}

	virtual int get_N_from_topology();
	virtual void read_topology(int N_from_conf, int *N_stars, BaseParticle<number> **particles);
	virtual void check_input_sanity(BaseParticle<number> **particles, int N);

	virtual void generate_random_configuration(BaseParticle<number> **particles, int N, number box_side);

	void set_only_intra(bool value) { _only_intra = value; }
};

#endif /* TSP_INTERACTION_H */
