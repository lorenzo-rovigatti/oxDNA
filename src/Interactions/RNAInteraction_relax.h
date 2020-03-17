/*
 * RNAInteraction_relax.h
 *
 *  Created on: Nov 28, 2014
 *      Author: Ben Snodin
 */

#ifndef RNA_INTERACTION_RELAX_H
#define RNA_INTERACTION_RELAX_H

#include "BaseInteraction.h"
#include "RNAInteraction.h"

/**
 * @brief Modified version of RNAInteraction which modifies the bonded backbone-backbone potential so it does not diverge
 *
 * Replaces the bonded backbone-backbone FENE potential with a harmonic potential. This is to allow very stressed initial
 * configurations, which might otherwise cause the simulation to fail, to relax to a sensible structure
 *
 * This interaction takes 3 compulsory arguments:
 *
 * This interaction is selected with
 * interaction_type = RNA_relax
 *
 @verbatim
 relax_type = <string> (Possible values: constant_force, harmonic_force; Relaxation algorithm used)
 relax_strength = <float> (Force constant for the replacement of the FENE potential)
 @endverbatim
 */

class RNAInteraction_relax: public RNAInteraction {
protected:
	inline virtual number _backbone(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	int _backbone_type;
	float _backbone_k;

	int _constant_force;
	int _harmonic_force;

public:
	RNAInteraction_relax();
	virtual ~RNAInteraction_relax();

	void check_input_sanity(std::vector<BaseParticle *> &particles);
	void get_settings(input_file &inp);
};

#endif /*RNA_INTERACTION_RELAX_H*/
