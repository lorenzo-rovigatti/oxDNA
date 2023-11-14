/*
 * DNANucleotide.cpp
 *
 *  Created on: 04/mar/2013
 *      Author: lorenzo
 */

#include "DNANucleotide.h"

LR_vector const DNANucleotide::principal_axis(1, 0, 0);
LR_vector const DNANucleotide::stack_axis(0, 0, 1);
LR_vector const DNANucleotide::third_axis(0, 1, 0);

DNANucleotide::DNANucleotide(bool grooving) :
				BaseParticle() {
	// _stack_axis definition has been changed for the major-minor grooving addition to the code and third_axis definition added;
	// _stack_axis is not used anywhere in the code apart from for major-minor grooving so the change should not affect anything
	// unfortunately the third axis does not in fact point along the z-axis which might be more intuitive.
	// -- Ben 1/8/13
	int_centers.resize(3);
	_grooving = grooving;
}

DNANucleotide::~DNANucleotide() {

}

bool DNANucleotide::is_bonded(BaseParticle *q) {
	return (n3 == q || n5 == q);
}

void DNANucleotide::set_positions() {
	if(_grooving) {
		int_centers[BACK] = orientation * principal_axis * POS_MM_BACK1 + orientation * third_axis * POS_MM_BACK2;
		int_centers[STACK] = orientation * principal_axis * POS_STACK;
	}
	else {
		int_centers[BACK] = orientation * principal_axis * POS_BACK;
		int_centers[STACK] = int_centers[BACK] * (POS_STACK / POS_BACK);
	}
	int_centers[BASE] = int_centers[STACK] * (POS_BASE / POS_STACK);
}
