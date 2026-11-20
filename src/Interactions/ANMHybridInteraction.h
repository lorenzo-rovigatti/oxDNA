/*
 * ANMHybridInteraction.h
 *
 *  Created Nov 2026
 *  Author: erik (from DNANMInteraction and RNANMInteraction by jonah and subhajit)
 *
 * Proteins described by the anisotropic network model (ANM) together with the nucleic acids of a nucleotide-level
 * interaction (NucleicInteraction = DNA2Interaction, RNA2Interaction or DRHInteraction).
 *
 * Every pair is sent to one of three places:
 * - nucleotide-nucleotide pairs go to NucleicInteraction, unchanged;
 * - amino acid-amino acid pairs go to the protein network model (ANMModel): springs, excluded volume and, in the
 *   "torsion" (ANMT) models, bending and torsion between consecutive amino acids;
 * - amino acid-nucleotide pairs get the ANMModel excluded volume between the amino acid and the backbone and base
 *   sites of the nucleotide.
 *
 * Subclasses only choose the nucleic acid type and how nucleotides are allocated.
 */

#ifndef ANMHYBRIDINTERACTION_H_
#define ANMHYBRIDINTERACTION_H_

#include "ANMModel.h"
#include "../Particles/BaseParticle.h"
#include "../Utilities/ANMTopology.h"

#include <string>
#include <vector>

template<class NucleicInteraction>
class ANMHybridInteraction: public NucleicInteraction {
protected:
	ANMModel _anm;
	std::string _nucleic_type;
	ANMTopology::Topology _anm_topology;
	/// true for the amino acids, indexed by particle index
	std::vector<bool> _is_amino_acid;

	/// allocates the nucleotide with the given index ('D' for DNA, 'R' for RNA)
	virtual BaseParticle *_new_nucleotide(int index, char kind) = 0;

	bool _is_protein(BaseParticle *p) const {
		return _is_amino_acid[p->index];
	}

	number _term_spring(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	number _term_protein_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	number _term_protein_nucleotide_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	number _term_angular(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	number _protein_nucleotide_excluded_volume(BaseParticle *p, BaseParticle *q, bool update_forces);

public:
	/// ids of the protein terms in the interaction map, after the ones of the nucleic acid interactions
	enum {
		SPRING = 8,
		PRO_EXC_VOL = 9,
		PRO_NUC_EXC_VOL = 10,
		PRO_ANG_POT = 11
	};

	/**
	 * @param torsion use the ANMT protein model (bending and torsion between consecutive amino acids)
	 * @param nucleic_type "DNA", "RNA" or "NA" (both, new-style topologies only)
	 */
	ANMHybridInteraction(bool torsion, const std::string &nucleic_type);
	virtual ~ANMHybridInteraction();

	void get_settings(input_file &inp) override;
	void init() override;

	void allocate_particles(std::vector<BaseParticle *> &particles) override;
	void read_topology(int *N_strands, std::vector<BaseParticle *> &particles) override;
	void check_input_sanity(std::vector<BaseParticle *> &particles) override;

	/// reads the topology counts without building the particles (used by the CUDA interactions, which are different
	/// objects from the one whose read_topology is called by the backend)
	void read_topology_counts();

	const ANMModel &anm() const {
		return _anm;
	}

	const ANMTopology::Topology &anm_topology() const {
		return _anm_topology;
	}

	number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false) override;
	number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false) override;
	number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false) override;
	number pair_interaction_term(int name, BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false) override;
};

#endif /* ANMHYBRIDINTERACTION_H_ */
