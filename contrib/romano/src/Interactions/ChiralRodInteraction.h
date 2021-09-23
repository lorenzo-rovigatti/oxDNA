/*
 * ChiralRodInteraction.h
 *
 *  Created on: 8/Oct/2013
 *      Author: Flavio
 */

#ifndef CHIRALRODINTERACTION_H_
#define CHIRALRODINTERACTION_H_

#include "Interactions/BaseInteraction.h"
#include "Interactions/InteractionUtils.h"
#include "Particles/SpheroCylinder.h"

/**
 * @brief Class to simulate chiral hard rods
 *
 *
 * The diameter of the cylinders is assumed to be 1.
 * The length is the length of the cylinder;
 * chiral_delta is the extent of the interaction;
 *
 * interaction_type = ChiralRod
 *
 * Input options:
 *
 * @verbatim
length = <float> (lenght of the cylinders)
chiral_delta = <float> (extent of the chiral interaction)
chiral_alpha = <float> (equilibrim angle)
chiral_interval = <float> (the interval for the chiral interaction)@endverbatim
 */

template <typename number>
class ChiralRodInteraction: public BaseInteraction<number, ChiralRodInteraction<number> > {
protected:
	/// length of the line segment
	number _length;

	number _chiral_delta;

	number _chiral_alpha;

	number _chiral_interval;

	number _chiral_max;

	number _chiral_min;

	number _epsilon;

	bool _attractive;

	number _att_range, _extra_range;

	number _min;

	/// potential
	inline number _chiral_pot (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);

public:
	enum {
		INNER_CYLINDER_OVERLAP = 0,
		OUTER_CYLINDER_OVERLAP = 1,
		ANGULAR_MODULATION = 2,
	};

	ChiralRodInteraction();
	virtual ~ChiralRodInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(BaseParticle<number> **particles, int N);

	virtual number pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_term(int name, BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, compute_r, update_forces);
	}

	virtual void check_input_sanity(BaseParticle<number> **particles, int N);

	bool generate_random_configuration_overlap (BaseParticle<number> * p, BaseParticle<number> *q);
};

extern "C" BaseInteraction<float> * make_float() { return new ChiralRodInteraction<float> (); }
extern "C" BaseInteraction<double> * make_double() { return new ChiralRodInteraction<double> (); }

#endif /* CHIRALRODINTERACTION_H_ */
