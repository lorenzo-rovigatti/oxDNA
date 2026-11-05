/*
 * DNANMInteraction.cpp
 *
 *  Created on: Feb 2026
 *  Updated Nov 2026
 *  Author: jonah, subhajit, erik
 */

#include "DNANMInteraction.h"

#include "../Particles/DNANucleotide.h"

DNANMInteraction::DNANMInteraction(bool torsion) :
				ANMHybridInteraction<DNA2Interaction>(torsion, "DNA") {

}

DNANMInteraction::~DNANMInteraction() = default;

BaseParticle *DNANMInteraction::_new_nucleotide(int index, char kind) {
	return new DNANucleotide(this->_grooving);
}
