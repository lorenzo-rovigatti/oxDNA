/*
 * DNAInteraction_relax2.h
 *
 *  Created on: Feb 25, 2016
 *      Author: Flavio
 */

#ifndef DNA_INTERACTION_RELAX2_H
#define DNA_INTERACTION_RELAX2_H

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
 * interaction_type = DNA_relax2
 *
@verbatim
relax_type = <string> (Possible values: constant_force, harmonic_force; Relaxation algorithm used)
relax_strength = <float> (Force constant for the replacement of the FENE potential)
@endverbatim
 */
template <typename number>
class DNAInteraction_relax2 : public DNAInteraction<number> {
protected:
	inline virtual number _backbone(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	number _fmax;

public:
	DNAInteraction_relax2();
	virtual ~DNAInteraction_relax2();

	void check_input_sanity(BaseParticle<number> **particles, int N);
	void get_settings(input_file &inp);
};

#endif /*DNA_INTERACTION_RELAX2_H*/
