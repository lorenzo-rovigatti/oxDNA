/*
 * JordanInteraction.h
 *
 *  Created on: 19/feb/2016
 *      Author: flavio
 */

#ifndef JORDANINTERACTION_H_
#define JORDANINTERACTION_H_

#include "BaseInteraction.h"

/**
 * @brief Manages the interaction between simple jordan particles (as described in http://jcp.aip.org/resource/1/jcpsa6/v131/i1/p014504_s1)
 *
 * This interaction is selected with
 * interaction_type = jordan
 *
 * @verbatim
JORDAN_N = <int> (number of patches)
[JORDAN_N_B = <int> (number of patches on species B)]
[JORDAN_alpha = <float> (width of patches, defaults to 0.12)]
@endverbatim
 */
template <typename number>
class JordanInteraction: public BaseInteraction<number, JordanInteraction<number> > {
protected:
	/// sigma of the gaussian modulation
	number _s;

	/// angle below the equator
	number _phi;

	/// "sensitivity" of the LJ contribution
	int _m;

	/// interaction
	number _jordan_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);


public:
	enum {
		JORDAN = 0
	};

	JordanInteraction();
	virtual ~JordanInteraction();

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

#endif /* JORDANINTERACTION_H_ */
