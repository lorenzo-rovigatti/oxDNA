
#ifndef POLYMER_INTERACTION_H
#define POLYMER_INTERACTION_H

#include <vector>
#include "Interactions/BaseInteraction.h"
#include "Particles/TSPParticle.h"
#include "Lists/Cells.h"

#include <numeric>

#define MAX_INSERTION_TRIES 1000000

template <typename number>
struct ChainDetails {
	int N() {
		return std::accumulate(block_lengths.begin(), block_lengths.end(), 0);
	}

	std::vector<int> block_lengths;
};

template <typename number>
class PolymerInteraction : public BaseInteraction<number, PolymerInteraction<number> > {
protected:
	number _rfene, _sqr_rfene;
	number _rfene_anchor, _sqr_rfene_anchor;

	number _Polymer_sqr_rep_rcut;
	number _Polymer_lambda;

	std::vector<ChainDetails<number>> _chains;

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
	PolymerInteraction();
	virtual ~PolymerInteraction();

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

	virtual void generate_random_configuration(BaseParticle<number> **particles, int N);
};

#endif /* POLYMER_INTERACTION_H */
