/*
 * ChargedCylInteraction.h
 *
 *  Created on: 10/Feb/2022
 *      Author: Flavio
 */

#ifndef CHARGEDCYLINTERACTION_H_
#define CHARGEDCYLINTERACTION_H_

#include "Interactions/BaseInteraction.h"
#include "Interactions/InteractionUtils.h"
#include "Particles/SpheroCylinder.h"

/**
 * @brief Class to simulate chiral hard rods
 *
 *
 * The diameter of the cylinders is assumed to be 1.
 * The length is the length of the cylinder;
 * lb is the characteristic length of the interaction
 *
 * interaction_type = ChiralRod
 *
 * Input options:
 *
 * @verbatim
length = <float> (lenght of the cylinders)
lb = <float> (extent of the chiral interaction)
 */

class ChargedCylInteraction: public BaseInteraction {
protected:
	/// length of the line segment
	number _length;

	/// characteristic length of the interaction
	number _lb;

	/// cutoff used for the charges
	number _charge_cutoff;

	/// distance between charges
	number _dc;

	/// number of charges
	int _n_charges;

	number _cylinder_overlap(BaseParticle *p, BaseParticle *q, bool _compute_r = true, bool update_forces=false);
	number _charge_interaction(BaseParticle *p, BaseParticle *q, bool _compute_r = true, bool update_forces=false);

public:
	enum {
		OVERLAP = 0,
		CHARGES = 1,
	};

	ChargedCylInteraction();
	virtual ~ChargedCylInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool _compute_r = true, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool _compute_r = true, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool _compute_r = true, bool update_forces=false);

    virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);

	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);

	//bool generate_random_configuration_overlap (BaseParticle * p, BaseParticle *q);
};

extern "C" ChargedCylInteraction * make_ChargedCylInteraction(); 

#endif /* CHARGEDCYLINTERACTION_H_ */
