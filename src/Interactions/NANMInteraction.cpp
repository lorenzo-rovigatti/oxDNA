/*
 * NANMInteraction.cpp
 *
 *  Created Nov 2026
 *  Author: erik
 */

#include "NANMInteraction.h"

#include "../Particles/DNANucleotide.h"
#include "../Particles/RNANucleotide.h"

NANMInteraction::NANMInteraction(bool torsion) :
				ANMHybridInteraction<DRHInteraction>(torsion, "NA") {

}

NANMInteraction::~NANMInteraction() = default;

void NANMInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
	// DRHInteraction uses _acid_type to tell DNA from RNA nucleotides. Amino acids are never passed to it
	_acid_type.assign(particles.size(), 'P');
	ANMHybridInteraction<DRHInteraction>::allocate_particles(particles);
}

BaseParticle *NANMInteraction::_new_nucleotide(int index, char kind) {
	_acid_type[index] = kind;
	if(kind == 'D') {
		return new DNANucleotide(DNA2Interaction::_grooving);
	}
	return new RNANucleotide();
}
