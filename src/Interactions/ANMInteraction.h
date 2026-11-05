/*
 * ANMInteraction.h
 *
 *  Created on: feb-2026
 *  Author: jonah, subhajit; restructured by erik Nov 2026
 */

#ifndef ANM_INTERACTION_H
#define ANM_INTERACTION_H

#include "BaseInteraction.h"
#include "ANMModel.h"
#include "../Utilities/ANMTopology.h"

/**
 * @brief Proteins described by the anisotropic network model (ANM), without nucleic acids.
 *
 * This interaction is selected with
 * interaction_type = AC (ANM) or ACT (ANMT: ANM plus bending and torsion between consecutive amino acids)
 *
 * Input options:
 *
 * @verbatim
 parfile = string (parameter file with the protein springs)
 massfile = string (set massfile for simulations w/ different massed items)
 bending_k = float (ACT only: global bending stiffness. If bending_k and torsion_k are not set they are read from the parameter file)
 torsion_k = float (ACT only: global torsion stiffness)
 topology = new-style topology (every strand is type=peptide) or the legacy protein format: <N> <N_strands>
 @endverbatim
 */
class ANMInteraction: public BaseInteraction {
protected:
	ANMModel _anm;

	number _spring(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	number _exc_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	number _angular(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

public:
	enum {
		SPRING_POTENTIAL = 0,
		EXC_VOL = 1,
		ANG_POT = 2
	};

	/// @param torsion use the ANMT protein model (ACT)
	explicit ANMInteraction(bool torsion = false);
	virtual ~ANMInteraction();

	void get_settings(input_file &inp) override;
	void init() override;

	void allocate_particles(std::vector<BaseParticle *> &particles) override;
	void read_topology(int *N_strands, std::vector<BaseParticle *> &particles) override;

	const ANMModel &anm() const {
		return _anm;
	}

	number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false) override;
	number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false) override;
	number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false) override;

	void check_input_sanity(std::vector<BaseParticle *> &particles) override;
};

#endif //ANM_INTERACTION_H
