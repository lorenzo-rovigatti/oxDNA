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
	if(_grooving == 2) {
		// QUESTION: since the type is going to be constant throughout the simulation, we probably
		// want to do this just once
		if(btype == 4) {
			_POS_MM_BACK1 = POS_MM_BACK1;
			_POS_MM_BACK2 = POS_MM_BACK2;
			_POS_STACK = POS_STACK;
			_POS_BASE = POS_BASE;
		}
		else {
			if(type == N_A) {
				_POS_MM_BACK1 = POS_MM_BACK1_A;
				_POS_MM_BACK2 = POS_MM_BACK2_A;
				_POS_STACK = POS_STACK_A;
				_POS_BASE = POS_BASE_A;
			}
			else if(type == N_G) {
				_POS_MM_BACK1 = POS_MM_BACK1_G;
				_POS_MM_BACK2 = POS_MM_BACK2_G;
				_POS_STACK = POS_STACK_G;
				_POS_BASE = POS_BASE_G;
			}
			else if(type == N_C) {
				_POS_MM_BACK1 = POS_MM_BACK1_C;
				_POS_MM_BACK2 = POS_MM_BACK2_C;
				_POS_STACK = POS_STACK_C;
				_POS_BASE = POS_BASE_C;
			}
			else if(type == N_T) {
				_POS_MM_BACK1 = POS_MM_BACK1_T;
				_POS_MM_BACK2 = POS_MM_BACK2_T;
				_POS_STACK = POS_STACK_T;
				_POS_BASE = POS_BASE_T;
			}
		}
		int_centers[BACK] = orientation * principal_axis * _POS_MM_BACK1 + orientation * third_axis * _POS_MM_BACK2;
		int_centers[STACK] = orientation * principal_axis * _POS_STACK;
		int_centers[BASE] = int_centers[STACK] * (_POS_BASE / _POS_STACK);

	}
	else if(_grooving == 1) {
		int_centers[BACK] = orientation * principal_axis * POS_MM_BACK1 + orientation * third_axis * POS_MM_BACK2;
		int_centers[STACK] = orientation * principal_axis * POS_STACK;
	    int_centers[BASE] = int_centers[STACK] * (POS_BASE / POS_STACK);

	}
	else if(_grooving == 0) {
		int_centers[BACK] = orientation * principal_axis * POS_BACK;
		int_centers[STACK] = int_centers[BACK] * (POS_STACK / POS_BACK);
        int_centers[BASE] = int_centers[STACK] * (POS_BASE / POS_STACK);
	}
}
