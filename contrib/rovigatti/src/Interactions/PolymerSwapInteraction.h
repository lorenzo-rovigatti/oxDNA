#ifndef POLYMERSWAP_INTERACTION_H
#define POLYMERSWAP_INTERACTION_H

#include <vector>
#include "Interactions/BaseInteraction.h"

/**
 * @brief Handles interactions in microgel systems.
 *
 * This interaction is selected with
 * interaction_type = PolymerSwapInteraction
 *
 * Input options:
 *
 * @verbatim

 @endverbatim
 */
class PolymerSwapInteraction: public BaseInteraction<PolymerSwapInteraction> {
protected:
	number _rfene, _sqr_rfene;

	number _PS_sqr_rep_rcut;
	number _PS_alpha, _PS_beta, _PS_gamma;
	int _PS_n;

	std::string _bond_filename;

	number _fene(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces);
	number _nonbonded(BaseParticle *p, BaseParticle *q, LR_vector *r, bool update_forces);

public:
	enum {
		BONDED = 0, NONBONDED = 1
	};
	PolymerSwapInteraction();
	virtual ~PolymerSwapInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false);
	virtual number pair_interaction_term(int name, BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, r, update_forces);
	}

	virtual void read_topology(int *N_stars, std::vector<BaseParticle *> &particles);
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
};

extern "C" PolymerSwapInteraction *make_PolymerSwapInteraction();

#endif /* POLYMERSWAP_INTERACTION_H */
