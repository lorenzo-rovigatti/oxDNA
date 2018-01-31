
#ifndef MG_INTERACTION_H
#define MG_INTERACTION_H

#include <vector>
#include "Interactions/BaseInteraction.h"

#define MAX_INSERTION_TRIES 1000000

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
template <typename number>
class MGInteraction : public BaseInteraction<number, MGInteraction<number> > {
protected:
	number _rfene, _sqr_rfene;

	number _MG_sqr_rep_rcut;
	number _MG_alpha, _MG_beta, _MG_gamma;
	int _MG_n;

	string _bond_filename;

	number _fene(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	number _nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);

public:
	enum {
		BONDED = 0,
		NONBONDED = 1
	};
	MGInteraction();
	virtual ~MGInteraction();

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
};

extern "C" MGInteraction<float> *make_MGInteraction_float();
extern "C" MGInteraction<double> *make_MGInteraction_double();

#endif /* MG_INTERACTION_H */
