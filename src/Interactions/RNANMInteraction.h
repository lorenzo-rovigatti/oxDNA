/**
 * @brief Hybrid RNA/ANM Model
 *
 * Jonah, Subhajit Feb. 2026; restructured by Erik Nov. 2026 (see ANMHybridInteraction)
 *
 * To use the oxRNA2 model with ANM proteins, set
 *
 * interaction_type = RNANM (or RNACT to add bending and torsion along the protein backbone)
 *
 * in the input file
 *
 * Input options:
 *
 * @verbatim
 parfile = string (parameter file with the protein springs, set to none for RNA-only simulations; mandatory for RNACT)
 massfile = string (set massfile for simulations w/ different massed items)
 bending_k = float (RNACT only: global bending stiffness. If bending_k and torsion_k are not set they are read from the parameter file)
 torsion_k = float (RNACT only: global torsion stiffness)
 topology = new-style topology (proteins are type=peptide strands) or the legacy mixed format:
            <N> <N_strands> <nrna> <npro> <nrnas>
 @endverbatim
 */

#ifndef RNANM_INTERACTION_H
#define RNANM_INTERACTION_H

#include "RNAInteraction2.h"
#include "ANMHybridInteraction.h"

class RNANMInteraction: public ANMHybridInteraction<RNA2Interaction> {
protected:
	BaseParticle *_new_nucleotide(int index, char kind) override;

public:
	/// @param torsion use the ANMT protein model (RNACT)
	explicit RNANMInteraction(bool torsion);
	virtual ~RNANMInteraction();
};

#endif /* RNANM_INTERACTION_H */
