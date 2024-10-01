/*
 * RNANucleotide.cpp
 *
 *  Created on: 04/dec/2013
 *      Author: petr
 */

#include "RNANucleotide.h"

RNANucleotide::RNANucleotide() :
				BaseParticle(),
				_principal_axis(1, 0, 0),
				_stack_axis(0, 0, 1),
				_third_axis(0, 1, 0) {
	int_centers.resize(7);
}

RNANucleotide::~RNANucleotide() {

}

bool RNANucleotide::is_bonded(BaseParticle *q) {
	return (n3 == q || n5 == q);
}

void RNANucleotide::set_positions() {
	_set_back_position();
	_set_stack_position();
	_set_base_position();
}

Model *RNANucleotide::model = 0;
