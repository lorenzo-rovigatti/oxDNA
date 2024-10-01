#ifndef MG_INTERACTION_H
#define MG_INTERACTION_H

#include <vector>
#include "Interactions/BaseInteraction.h"

/**
 * @brief Handles interactions in microgel systems.
 *
 * This interaction is selected with
 * interaction_type = MGInteraction
 *
 * Input options:
 *
 * @verbatim

 @endverbatim
 */
class MGInteraction: public BaseInteraction {
protected:
	number _rfene, _sqr_rfene;

	number _MG_sqr_rep_rcut;
	number _MG_alpha, _MG_beta, _MG_gamma;
	int _MG_n;

	std::string _bond_filename;

	number _fene(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	number _nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

public:
	enum {
		BONDED = 0, NONBONDED = 1
	};
	MGInteraction();
	virtual ~MGInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	virtual void read_topology(int *N_stars, std::vector<BaseParticle *> &particles);
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
};

extern "C" MGInteraction *make_MGInteraction();

#endif /* MG_INTERACTION_H */
