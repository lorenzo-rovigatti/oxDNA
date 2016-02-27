/*
 * PatchyParticle.cpp
 *
 *  Created on: 15/mar/2013
 *      Author: lorenzo
 */

#include "PatchyParticle.h"
#include "../Utilities/oxDNAException.h"

#define HALF_ISQRT3 0.28867513459481292f

template<typename number>
PatchyParticle<number>::PatchyParticle(int N_patches, int nt, number sigma) : BaseParticle<number>(), _sigma(sigma) {
	this->type = nt;
	this->N_int_centers = N_patches;
	this->int_centers = new LR_vector<number>[N_patches];
	_base_patches = new LR_vector<number>[N_patches];

	_set_base_patches();
}

template<typename number>
PatchyParticle<number>::~PatchyParticle() {
	delete[] _base_patches;
}

template<typename number>
void PatchyParticle<number>::_set_base_patches() {
	switch(this->N_int_centers) {
	case 2: {
		_base_patches[0] = LR_vector<number>(0, 1, 0);
		_base_patches[1] = LR_vector<number>(0, -1, 0);

		break;
	}
	case 3: {
		number cos30 = cos(M_PI / 6.);
		number sin30 = sin(M_PI / 6.);

		_base_patches[0] = LR_vector<number>(0, 1, 0);
		_base_patches[1] = LR_vector<number>(cos30, -sin30, 0);
		_base_patches[2] = LR_vector<number>(-cos30, -sin30, 0);

		break;
	}
	case 4: {
		_base_patches[0] = LR_vector<number>(-HALF_ISQRT3, -HALF_ISQRT3,  HALF_ISQRT3);
		_base_patches[1] = LR_vector<number>( HALF_ISQRT3, -HALF_ISQRT3, -HALF_ISQRT3);
		_base_patches[2] = LR_vector<number>( HALF_ISQRT3,  HALF_ISQRT3,  HALF_ISQRT3);
		_base_patches[3] = LR_vector<number>(-HALF_ISQRT3,  HALF_ISQRT3, -HALF_ISQRT3);
		break;
	}
	default:
		throw oxDNAException("Unsupported number of patches %d\n", this->N_int_centers);
	}

	for(int i = 0; i < this->N_int_centers; i++) {
		_base_patches[i].normalize();
		_base_patches[i] *= 0.5;
	}
}

template<typename number>
void PatchyParticle<number>::set_positions() {
	for(int i = 0; i < this->N_int_centers; i++) this->int_centers[i] = (this->orientation*_base_patches[i])*_sigma;
}

template class PatchyParticle<float>;
template class PatchyParticle<double>;
