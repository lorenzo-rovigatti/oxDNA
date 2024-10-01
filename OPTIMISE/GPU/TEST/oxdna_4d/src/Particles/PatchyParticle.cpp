/*
 * PatchyParticle.cpp
 *
 *  Created on: 15/mar/2013
 *      Author: lorenzo
 */

#include "PatchyParticle.h"
#include "../Utilities/oxDNAException.h"

#define HALF_ISQRT3 0.28867513459481292f

PatchyParticle::PatchyParticle(std::vector<LR_vector> base_patches, int nt, number sigma) :
				BaseParticle(),
				_sigma(sigma),
				_base_patches(base_patches),
				_N_patches(base_patches.size()) {
	type = btype = nt;
	int_centers.resize(base_patches.size());

	for(uint i = 0; i < N_int_centers(); i++) {
		_base_patches[i].normalize();
		_base_patches[i] *= 0.5;
	}
}

PatchyParticle::PatchyParticle(int N_patches, int nt, number sigma) :
				BaseParticle(),
				_sigma(sigma),
				_N_patches(N_patches) {
	type = btype = nt;
	int_centers.resize(N_patches);
	_base_patches.resize(N_patches);

	switch(N_int_centers()) {
	case 0:
		break;
	case 1: {
		_base_patches[0] = LR_vector(1, 0, 0);
		break;
	}
	case 2: {
		_base_patches[0] = LR_vector(0, 1, 0);
		_base_patches[1] = LR_vector(0, -1, 0);

		break;
	}
	case 3: {
		number cos30 = cos(M_PI / 6.);
		number sin30 = sin(M_PI / 6.);

		_base_patches[0] = LR_vector(0, 1, 0);
		_base_patches[1] = LR_vector(cos30, -sin30, 0);
		_base_patches[2] = LR_vector(-cos30, -sin30, 0);

		break;
	}
	case 4: {
		_base_patches[0] = LR_vector(-HALF_ISQRT3, -HALF_ISQRT3, HALF_ISQRT3);
		_base_patches[1] = LR_vector( HALF_ISQRT3, -HALF_ISQRT3, -HALF_ISQRT3);
		_base_patches[2] = LR_vector( HALF_ISQRT3, HALF_ISQRT3, HALF_ISQRT3);
		_base_patches[3] = LR_vector(-HALF_ISQRT3, HALF_ISQRT3, -HALF_ISQRT3);
		break;
	}
	case 5: {
		_base_patches[0] = LR_vector(0.135000, -0.657372, -0.741375);
		_base_patches[1] = LR_vector(0.259200, 0.957408, -0.127224);
		_base_patches[2] = LR_vector(-0.394215, -0.300066, 0.868651);
		_base_patches[3] = LR_vector(-0.916202, 0.202077, -0.346033);
		_base_patches[4] = LR_vector(0.916225, -0.202059, 0.345982);
		break;
	}
	default:
		throw oxDNAException("Unsupported number of patches %d\n", N_int_centers());
	}

	for(uint i = 0; i < N_int_centers(); i++) {
		_base_patches[i].normalize();
		_base_patches[i] *= 0.5;
	}
}

PatchyParticle::~PatchyParticle() {

}

void PatchyParticle::set_positions() {
	for(uint i = 0; i < N_int_centers(); i++) {
		int_centers[i] = (orientation * _base_patches[i]) * _sigma;
	}
}
