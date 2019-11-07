#ifndef POLYMER_INTERACTION_H
#define POLYMER_INTERACTION_H

#include <vector>
#include "Interactions/BaseInteraction.h"
#include "Particles/TSPParticle.h"
#include "Lists/Cells.h"

#include <numeric>

struct ChainDetails {
	int N() {
		return std::accumulate(block_lengths.begin(), block_lengths.end(), 0);
	}

	std::vector<int> block_lengths;
};

class PolymerInteraction: public BaseInteraction<PolymerInteraction> {
protected:
	number _rfene, _sqr_rfene;

	number _Polymer_sqr_rep_rcut;
	number _Polymer_lambda;
	number _Polymer_lambda_E_cut;

	std::vector<ChainDetails> _chains;

	number _fene(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces);
	number _nonbonded(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces);

public:
	enum {
		BONDED = 0, NONBONDED = 1
	};
	PolymerInteraction();
	virtual ~PolymerInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(BaseParticle **particles, int N);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false);
	virtual number pair_interaction_term(int name, BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, r, update_forces);
	}

	virtual int get_N_from_topology();
	virtual void read_topology(int N_from_conf, int *N_stars, BaseParticle **particles);
	virtual void check_input_sanity(BaseParticle **particles, int N);

	virtual void generate_random_configuration(BaseParticle **particles, int N);
};

extern "C" PolymerInteraction *make_PolymerInteraction();

#endif /* POLYMER_INTERACTION_H */
