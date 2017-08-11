/*
 * DNA2Interaction_relax.h
 *
 *  Created on: Mar 13, 2017
 *      Author: Petr Sulc
 */

#ifndef DNA2_INTERACTION_RELAX_H
#define DNA2_INTERACTION_RELAX_H

#include "BaseInteraction.h"
#include "DNA2Interaction.h"

// HARMONIC BONDED BACKBONE-BACKBONE
#define HARMONIC2_R0 FENE_R0_OXDNA2

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
template <typename number>
class DNA2Interaction_relax : public DNA2Interaction<number> {
protected:
	inline virtual number _backbone(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	int _backbone_type;
	float _backbone_k;

	int _constant_force;
	int _harmonic_force;

public:
	DNA2Interaction_relax();
	virtual ~DNA2Interaction_relax();

	void check_input_sanity(BaseParticle<number> **particles, int N);
	void get_settings(input_file &inp);
};

#endif /*DNA2_INTERACTION_RELAX_H*/
