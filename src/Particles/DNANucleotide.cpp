/*
 * DNANucleotide.cpp
 *
 *  Created on: 04/mar/2013
 *      Author: lorenzo
 */

#include "DNANucleotide.h"
#include <iostream>

LR_vector const DNANucleotide::principal_axis(1, 0, 0);
LR_vector const DNANucleotide::stack_axis(0, 0, 1);
LR_vector const DNANucleotide::third_axis(0, 1, 0);

const std::array<number, 5> DNANucleotide::all_pos_mm_back1({
	POS_MM_BACK1,
	POS_MM_BACK1_A,
	POS_MM_BACK1_G,
	POS_MM_BACK1_C,
	POS_MM_BACK1_T
});
const std::array<number, 5> DNANucleotide::all_pos_mm_back2({
	POS_MM_BACK2,
	POS_MM_BACK2_A,
	POS_MM_BACK2_G,
	POS_MM_BACK2_C,
	POS_MM_BACK2_T
});
const std::array<number, 5> DNANucleotide::all_pos_stack({
	POS_STACK,
	POS_STACK_A,
	POS_STACK_G,
	POS_STACK_C,
	POS_STACK_T
});
const std::array<number, 5> DNANucleotide::all_pos_base({
	POS_BASE,
	POS_BASE_A,
	POS_BASE_G,
	POS_BASE_C,
	POS_BASE_T
});

DNANucleotide::DNANucleotide(int grooving) :
				BaseParticle()
{
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
	if(_grooving > 0) { // DNA2 and DNA3
		if(pos_mm_back1 == 0.) { // we do this just once at the beginning
			// for DNA2, always use idx=0
			// for DNA3, map btype==4 to idx=0, otherwise idx=1+type (recalling that N_A=0, N_G=1, N_C=2, N_T=3)
			int idx = (_grooving == 1 || btype == 4) ? 0 : 1 + type;
			pos_mm_back1 = all_pos_mm_back1[idx];
			pos_mm_back2 = all_pos_mm_back2[idx];
			pos_stack = all_pos_stack[idx];
			pos_base = all_pos_base[idx];
		}

		// branchless
		int_centers[BACK] = orientation * principal_axis * pos_mm_back1 + orientation * third_axis * pos_mm_back2;
		int_centers[STACK] = orientation * principal_axis * pos_stack;
		int_centers[BASE] = int_centers[STACK] * (pos_base / pos_stack);
	}
	else if(_grooving == 0) { // DNA1
		int_centers[BACK] = orientation * principal_axis * POS_BACK;
		int_centers[STACK] = int_centers[BACK] * (POS_STACK / POS_BACK);
        int_centers[BASE] = int_centers[STACK] * (POS_BASE / POS_STACK);
	}
}
