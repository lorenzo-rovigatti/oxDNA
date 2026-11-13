/**
 * @brief Hybrid DNA/RNA/ANM model: the oxNA (NA) interaction for DNA, RNA and DNA/RNA hybrids together with ANM proteins.
 *
 * Erik Nov. 2026 (see ANMHybridInteraction)
 *
 * To use it, set
 *
 * interaction_type = NANM (or NACT to add bending and torsion along the protein backbone)
 *
 * in the input file. As for the NA interaction, the topology must be in the new (5->3) format, with each strand
 * having type=DNA, type=RNA or type=peptide.
 *
 * Input options: those of the NA interaction, plus
 *
 * @verbatim
 parfile = string (parameter file with the protein springs, set to none for simulations without proteins; mandatory for NACT)
 massfile = string (set massfile for simulations w/ different massed items)
 bending_k = float (NACT only: global bending stiffness. If bending_k and torsion_k are not set they are read from the parameter file)
 torsion_k = float (NACT only: global torsion stiffness)
 @endverbatim
 */

#ifndef NANM_INTERACTION_H
#define NANM_INTERACTION_H

#include "DRHInteraction.h"
#include "ANMHybridInteraction.h"

class NANMInteraction: public ANMHybridInteraction<DRHInteraction> {
protected:
	BaseParticle *_new_nucleotide(int index, char kind) override;

public:
	/// @param torsion use the ANMT protein model (NACT)
	explicit NANMInteraction(bool torsion);
	virtual ~NANMInteraction();

	void allocate_particles(std::vector<BaseParticle *> &particles) override;
};

#endif /* NANM_INTERACTION_H */
