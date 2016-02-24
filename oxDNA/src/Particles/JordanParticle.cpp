/*
 * JordanParticle.cpp
 *
 *  Created on: 15/mar/2013
 *      Author: lorenzo
 */

#include "JordanParticle.h"
#include "../Utilities/oxDNAException.h"


template<typename number>
JordanParticle<number>::JordanParticle(number phi) : BaseParticle<number>() {
	int N_patches = 3;
	this->N_int_centers = N_patches;
	this->int_centers = new LR_vector<number>[N_patches];
	_base_patches = new LR_vector<number>[N_patches];

	_set_base_patches(phi);
}

template<typename number>
JordanParticle<number>::~JordanParticle() {
	delete[] _base_patches;
}

template<typename number>
void JordanParticle<number>::_set_base_patches(number phi) {
	switch(this->N_int_centers) {
		/*
		case 2: {
			_base_patches[0] = LR_vector<number>(0, 1, 0);
			_base_patches[1] = LR_vector<number>(0, -1, 0);
			break;
		}
		*/
		case 3: {
			_base_patches[0] = LR_vector<number>(cos(phi), 0, -sin(phi));
			_base_patches[1] = LR_vector<number>(cos(2.*M_PI/3.)*cos(phi), sin(2.*M_PI/3.)*cos(phi), -sin(phi));
			_base_patches[2] = LR_vector<number>(cos(4.*M_PI/4.)*cos(phi), sin(4.*M_PI/3.)*cos(phi), -sin(phi));

			break;
		}
		default:
			throw oxDNAException("Unsupported number of patches %d\n", this->N_int_centers);
	}

	for(int i = 0; i < this->N_int_centers; i++) {
		_base_patches[i].normalize();
		_base_patches[i] *= 0.5;
	}

	set_positions();
}

template<typename number>
void JordanParticle<number>::set_positions() {
	for(int i = 0; i < this->N_int_centers; i++) this->int_centers[i] = this->orientation * _base_patches[i];
}

template class JordanParticle<float>;
template class JordanParticle<double>;
