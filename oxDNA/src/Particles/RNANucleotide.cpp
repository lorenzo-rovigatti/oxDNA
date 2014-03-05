/*
 * RNANucleotide.cpp
 *
 *  Created on: 04/dec/2013
 *      Author: petr
 */

#include "RNANucleotide.h"

template<typename number>
RNANucleotide<number>::RNANucleotide() : BaseParticle<number>(), _principal_axis(1,0,0), _stack_axis(0,0,1), _third_axis(0,1,0) {
	//
	this->int_centers = new LR_vector<number>[7];
	this->N_int_centers = 7;
}

template<typename number>
RNANucleotide<number>::~RNANucleotide() {

}

template<typename number>
bool RNANucleotide<number>::is_bonded(BaseParticle<number> *q) {
	return (this->n3 == q || this->n5 == q);
}

template<typename number>
void RNANucleotide<number>::set_positions() {
	_set_back_position();
	_set_stack_position();
	_set_base_position();
}

template class RNANucleotide<double>;
template class RNANucleotide<float>;

template<class number> Model* RNANucleotide<number>::model = 0;
//template class Model* RNANucleotide<double>::model = 0;
