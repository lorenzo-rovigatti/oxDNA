/*
 * ChiralRodExplicit.h
 *
 *  Created on: 8/Oct/2013
 *      Author: Flavio
 */

#ifndef CHIRALRODEXPLICIT_H_
#define CHIRALRODEXPLICIT_H_

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
@verbatim
length = <float> (lenght of the cylinders)
chiral_delta = <float> (extent of the chiral interaction)
chiral_alpha = <float> (equilibrim angle)
chiral_interval = <float> (the interval for the chiral interaction)
sigma_solvent = <float> (radius of the solvent colloids)
@endverbatim
 */

template <typename number>
class ChiralRodExplicit: public BaseInteraction<number, ChiralRodExplicit<number> > {
protected:
	/// length of the line segment
	number _length;

	/// radial extent of the chiral interaction
	number _chiral_delta;

	/// preferred angle for chiral interactions
	number _chiral_alpha;

	/// angular acceptance of the chiral interaction
	number _chiral_interval;

	/// helper values computed in init: maximum and minumum angle for the chiral interaction
	number _chiral_max, _chiral_min;

	/// radius of the solvent hard spheres
	number _sigma_solv;

	/// number of chiral cylinders and of solvent molecules; read in get_settings
	int _N_rods, _N_solv;

	bool _togrow;

	/// interaction between chiral rods
	inline number _chiral_pot (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);

	/// interaction between rods and solvent
	inline number _rod_solv (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);

	/// interaction between solvent and solvent
	inline number _solv_solv (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);


public:
	enum {
		INNER_CYLINDER_OVERLAP = 0,
		OUTER_CYLINDER_OVERLAP = 1,
		ANGULAR_MODULATION = 2,
	};

	ChiralRodExplicit();
	virtual ~ChiralRodExplicit();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(BaseParticle<number> **particles, int N);

	virtual number pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_term(int name, BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, compute_r, update_forces);
	}

	virtual void read_topology(int *N_strands, BaseParticle<number> **particles);
	virtual void check_input_sanity(BaseParticle<number> **particles, int N);

	bool generate_random_configuration_overlap (BaseParticle<number> * p, BaseParticle<number> *q);
};

extern "C" BaseInteraction<float> * make_float() { return new ChiralRodExplicit<float> (); }
extern "C" BaseInteraction<double> * make_double() { return new ChiralRodExplicit<double> (); }

#endif /* CHIRALRODEXPLICIT_H_ */
