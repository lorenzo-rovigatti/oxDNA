/**
 * @brief Hybrid DNA/ANM Model
 *
 * Jonah, Subhajit Feb. 2026; restructured by Erik Nov. 2026 (see ANMHybridInteraction)
 *
 * To use the oxDNA2 model with ANM proteins, set
 *
 * interaction_type = DNANM (or DNACT to add bending and torsion along the protein backbone)
 *
 * in the input file
 *
 * Input options:
 *
 * @verbatim
 parfile = string (parameter file with the protein springs, set to none for DNA-only simulations; mandatory for DNACT)
 massfile = string (set massfile for simulations w/ different massed items)
 bending_k = float (DNACT only: global bending stiffness. If bending_k and torsion_k are not set they are read from the parameter file)
 torsion_k = float (DNACT only: global torsion stiffness)
 topology = new-style topology (proteins are type=peptide strands) or the legacy mixed format:
            <N> <N_strands> <ndna> <npro> <ndnas>
 @endverbatim
 */

#ifndef DNANM_INTERACTION_H
#define DNANM_INTERACTION_H

#include "DNA2Interaction.h"
#include "ANMHybridInteraction.h"

class DNANMInteraction: public ANMHybridInteraction<DNA2Interaction> {
protected:
	BaseParticle *_new_nucleotide(int index, char kind) override;

public:
	/// @param torsion use the ANMT protein model (DNACT)
	explicit DNANMInteraction(bool torsion);
	virtual ~DNANMInteraction();
};

#endif /* DNANM_INTERACTION_H */
