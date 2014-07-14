/*
 * DNANucleotide.cpp
 *
 *  Created on: 04/mar/2013
 *      Author: lorenzo
 */

#include "DNANucleotide.h"

template <typename number> LR_vector<number> const DNANucleotide<number>::principal_axis(1, 0, 0);
template <typename number> LR_vector<number> const DNANucleotide<number>::stack_axis(0, 0, 1);
template <typename number> LR_vector<number> const DNANucleotide<number>::third_axis(0, 1, 0);

template<typename number>
DNANucleotide<number>::DNANucleotide(bool grooving) : BaseParticle<number>() {
	// _stack_axis definition has been changed for the major-minor grooving addition to the code and third_axis definition added;
	// _stack_axis is not used anywhere in the code apart from for major-minor grooving so the change should not affect anything
	// unfortunately the third axis does not in fact point along the z-axis which might be more intuitive.
	// -- Ben 1/8/13
	this->int_centers = new LR_vector<number>[3];
	this->N_int_centers = 3;
	_grooving = grooving;
}

template<typename number>
DNANucleotide<number>::~DNANucleotide() {

}

template<typename number>
bool DNANucleotide<number>::is_bonded(BaseParticle<number> *q) {
	return (this->n3 == q || this->n5 == q);
}

template<typename number>
void DNANucleotide<number>::set_positions() {
	if (_grooving) {
		this->int_centers[BACK] = this->orientation*principal_axis*POS_MM_BACK1 + this->orientation*third_axis*POS_MM_BACK2;
		this->int_centers[STACK] = this->orientation*principal_axis*POS_STACK;
	}
	else{
		this->int_centers[BACK] = this->orientation*principal_axis*POS_BACK;
		this->int_centers[STACK] = this->int_centers[BACK]*(POS_STACK/POS_BACK);
	}
	this->int_centers[BASE] = this->int_centers[STACK]*(POS_BASE/POS_STACK);
}

template class DNANucleotide<double>;
template class DNANucleotide<float>;
