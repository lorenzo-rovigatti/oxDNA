/*
 * JordanParticle.cpp
 *
 *  Created on: 15/mar/2013
 *      Author: lorenzo
 */

#include "JordanParticle.h"
#include "../Utilities/oxDNAException.h"



JordanParticle::JordanParticle(int npatches, number phi, number int_k) : BaseParticle() {
	int N_patches = npatches;
	int_centers.resize(N_patches);
	_base_patches = new LR_vector[N_patches];
	_patch_rotations = new LR_matrix[N_patches];
	
	_set_base_patches(phi);
	_int_k = int_k;
}


JordanParticle::~JordanParticle() {
	delete[] _base_patches;
	delete[] _patch_rotations;
}


void JordanParticle::_set_base_patches(number phi) {
	switch(N_int_centers()) {
		/*
		case 2: {
			_base_patches[0] = LR_vector(0, 1, 0);
			_base_patches[1] = LR_vector(0, -1, 0);
			break;
		}
		*/
		case 3: {
			_base_patches[0] = LR_vector(cos(phi), 0, -sin(phi));
			_base_patches[1] = LR_vector(cos(2.*M_PI/3.)*cos(phi), sin(2.*M_PI/3.)*cos(phi), -sin(phi));
			_base_patches[2] = LR_vector(cos(4.*M_PI/4.)*cos(phi), sin(4.*M_PI/3.)*cos(phi), -sin(phi));

			break;
		}
		case 4: {
			_base_patches[0] = LR_vector(cos(phi), 0, -sin(phi));
			_base_patches[1] = LR_vector(cos(1.*M_PI/2.)*cos(phi), sin(1.*M_PI/2.)*cos(phi), -sin(phi));
			_base_patches[2] = LR_vector(cos(2.*M_PI/2.)*cos(phi), sin(2.*M_PI/2.)*cos(phi), -sin(phi));
			_base_patches[3] = LR_vector(cos(3.*M_PI/2.)*cos(phi), sin(3.*M_PI/2.)*cos(phi), -sin(phi));

			break;
		}
		default:
			throw oxDNAException("Unsupported number of patches %d\n", N_int_centers());
	}

	for(uint i = 0; i < N_int_centers(); i++) {
		_base_patches[i].normalize();
		_base_patches[i] *= 0.5;
		_patch_rotations[i] = LR_matrix (1., 0., 0.,  0., 1., 0.,  0., 0., 1.);
	}

	set_positions();
}


number JordanParticle::int_potential() {
	number energy = 0.f;
	
	for(uint i = 0; i < N_int_centers(); i ++) {
		// the trace of a rotation matrix == 1. + 2.cos(t)
		number arg = _patch_rotations[i].v1.x + _patch_rotations[i].v2.y + _patch_rotations[i].v3.z;
		arg = arg - 1.f;
		if (fabs(arg) > 1.f) arg = copysign(1.f, arg);
		number gamma = acos(arg);
		energy += _int_k * gamma * gamma;
	}
	
	return energy;
}


void JordanParticle::set_positions() {
	for(uint i = 0; i < N_int_centers(); i++) {
		this->int_centers[i] = this->orientation * (_patch_rotations[i] * _base_patches[i]);
	}
}


std::string JordanParticle::get_output_string () {
	std::string ret = Utils::sformat("%g %g %g ", this->pos.x, this->pos.y, this->pos.z);
	for (uint i = 0; i < N_int_centers(); i ++) {
		LR_vector patch = _patch_rotations[i] * _base_patches[i]; 
		ret += Utils::sformat("%g %g %g ", patch.x, patch.y, patch.z);
	}
	return ret;
}
