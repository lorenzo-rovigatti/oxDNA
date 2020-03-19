#ifndef POLYMERSWAP_INTERACTION_H
#define POLYMERSWAP_INTERACTION_H

#include "Interactions/BaseInteraction.h"

#include <vector>
#include <array>

/**
 * @brief Handles interactions in microgel systems.
 *
 * This interaction is selected with
 * interaction_type = PolymerSwapInteraction
 *
 * @verbatim

 @endverbatim
 */
class PolymerSwapInteraction: public BaseInteraction<PolymerSwapInteraction> {
protected:
	std::array<number, 3> _Kfene = {15., 15., 15.};
	std::array<number, 3> _rfene = {1.5, 1.0, 1.5};
	std::array<number, 3> _sqr_rfene = {2.25, 1.0, 2.25};
	std::array<number, 3> _WCA_sigma = {1.0, 1.0, 1.0};
	std::array<number, 3> _PS_sqr_rep_rcut;

	number _PS_alpha = 0.;
	number _PS_beta = 0.;
	number _PS_gamma = 0.;
	int _PS_n = 6;

	std::string _bond_filename;

	bool _only_links_in_bondfile = true;

	number _fene(BaseParticle *p, BaseParticle *q, bool update_forces);
	number _nonbonded(BaseParticle *p, BaseParticle *q, bool update_forces);

public:
	enum {
		BONDED = 0, NONBONDED = 1
	};
	enum {
		MONOMER = 0, STICKY = 1
	};
	PolymerSwapInteraction();
	virtual ~PolymerSwapInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_term(int name, BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false) {
		return _pair_interaction_term_wrapper(this, name, p, q, compute_r, update_forces);
	}

	virtual void read_topology(int *N_stars, std::vector<BaseParticle *> &particles);
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
};

extern "C" PolymerSwapInteraction *make_PolymerSwapInteraction();

#endif /* POLYMERSWAP_INTERACTION_H */
