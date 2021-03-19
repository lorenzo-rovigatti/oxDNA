/*
 * PatchyParticleDan.cpp
 *
 *  Created on: 01/feb/2016
 *      Author: lorenzo -> dan
 */

#include "PatchyParticleDan.h"
#include "../Utilities/oxDNAException.h"

PatchyParticleDan::PatchyParticleDan(int _N_patches, LR_vector * _inp_patch_vectors, LR_vector * _inp_ref_vectors, bool tor_flag) :
				BaseParticle() {

	_base_patch_vectors = NULL;
	_base_ref_vectors = NULL;
	_ref_vectors = NULL;

	//Create new arrays (one row for each patch)
	int_centers.resize(_N_patches);
	//printf("PP 3\n");
	_base_patch_vectors = new LR_vector[_N_patches];
	//printf("PP 4\n");
	_base_ref_vectors = new LR_vector[_N_patches];
	_ref_vectors = new LR_vector[_N_patches];
	//printf("PP 5\n");

	//Define _base_patch_vectors and _base_ref_vectors
	for(uint patch = 0; patch < this->N_int_centers(); patch++) {
		_base_patch_vectors[patch] = _inp_patch_vectors[patch];
		_base_ref_vectors[patch] = _inp_ref_vectors[patch];
		//printf("PP 6 patch %d _base_patch_vectors %f %f %f\n", patch, _base_patch_vectors[patch].x, _base_patch_vectors[patch].y, _base_patch_vectors[patch].z);

	}

	//printf("PP 7\n");

	/*for(int patch = 0; patch < this->N_int_centers(); patch++) {
	 printf("PP patch %d _base_patch_vectors %f %f %f, _base_ref_vectors %f %f %f\n", patch, _base_patch_vectors[patch].x, _base_patch_vectors[patch].y, _base_patch_vectors[patch].z, _base_ref_vectors[patch].x, _base_ref_vectors[patch].y, _base_ref_vectors[patch].z);
	 }*/

	//Set lengths of _base_patch_vectors (0.5) and _base_ref_vectors (1.0)
	for(int patch = 0; patch < _N_patches; patch++) {
		_base_patch_vectors[patch].normalize();
		_base_patch_vectors[patch] *= 0.5;
		//printf("PP 8 patch %d _base_patch_vectors %f %f %f\n", patch, _base_patch_vectors[patch].x, _base_patch_vectors[patch].y, _base_patch_vectors[patch].z);
		_base_ref_vectors[patch].normalize();

	}

	//printf("PP 9\n");
}

PatchyParticleDan::~PatchyParticleDan() {
	//printf("PP, ~PatchyParticleDan\n");

	/*Not sure if this is needed, and if it is, should it follow format below? 10/11/16
	 delete[] this->int_centers;*/
	if(_base_patch_vectors != NULL)
		delete[] _base_patch_vectors;
	if(_base_ref_vectors != NULL)
		delete[] _base_ref_vectors;
	if(_ref_vectors != NULL)
		delete[] _ref_vectors;

}

void PatchyParticleDan::set_positions() {
	//Set patch vectors and reference vectors based on particle orientation
	for(uint patch = 0; patch < this->N_int_centers(); patch++) {
		this->int_centers[patch] = this->orientation * _base_patch_vectors[patch];
		_ref_vectors[patch] = this->orientation * _base_ref_vectors[patch];
	}

	//printf("PP 10\n");
}

void PatchyParticleDan::copy_from(const BaseParticle &p) {
	//To inherit/extend copy_from from BaseParticle (see Lorenzo email)
	BaseParticle::copy_from(p);
	//To add the same functionality for ref_vectors

	const auto ppp = dynamic_cast<const PatchyParticleDan *>(&p);
	for(uint patch = 0; patch < this->N_int_centers(); patch++) {
		_ref_vectors[patch] = ppp->_ref_vectors[patch];
	}
}
