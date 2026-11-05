/*
 * RNANMInteraction.cpp
 *
 *  Created on: 2026
 *  Updated Nov 2026
 *  Author: jonah, subhajit, erik
 */

#include "RNANMInteraction.h"

#include "../Particles/RNANucleotide.h"

RNANMInteraction::RNANMInteraction(bool torsion) :
				ANMHybridInteraction<RNA2Interaction>(torsion, "RNA") {

}

RNANMInteraction::~RNANMInteraction() = default;

BaseParticle *RNANMInteraction::_new_nucleotide(int index, char kind) {
	return new RNANucleotide();
}
