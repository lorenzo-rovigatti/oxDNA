/*
 * DNAInteraction_relax.h
 *
 *  Created on: Mar 13, 2013
 *      Author: Ben Snodin
 */

#ifndef DNA_INTERACTION_RELAX_H
#define DNA_INTERACTION_RELAX_H

#include "BaseInteraction.h"
#include "DNAInteraction.h"

// HARMONIC BONDED BACKBONE-BACKBONE
#define HARMONIC_R0 FENE_R0_OXDNA

/**
 * @brief Modified version of DNAInteraction which modifies the bonded backbone-backbone potential so it does not diverge
 *
 * Replaces the bonded backbone-backbone FENE potential with a harmonic potential. This is to allow very stressed initial
 * configurations, which might otherwise cause the simulation to fail, to relax to a sensible structure
 *
 * This interaction takes 3 compulsory arguments:
 *
 * This interaction is selected with
 * interaction_type = DNA_relax
 *
 @verbatim
 relax_type = <string> (Possible values: constant_force, harmonic_force; Relaxation algorithm used)
 relax_strength = <float> (Force constant for the replacement of the FENE potential)
 @endverbatim
 */

class DNAInteraction_relax: public DNAInteraction {
protected:
	inline virtual number _backbone(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	int _backbone_type;
	float _backbone_k;

	int _constant_force;
	int _harmonic_force;

public:
	DNAInteraction_relax();
	virtual ~DNAInteraction_relax();

	void check_input_sanity(std::vector<BaseParticle *> &particles);
	void get_settings(input_file &inp);
};

#endif /*DNA_INTERACTION_RELAX_H*/
