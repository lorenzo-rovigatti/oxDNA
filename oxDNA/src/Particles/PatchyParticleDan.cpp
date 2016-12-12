/*
 * PatchyParticleDan.cpp
 *
 *  Created on: 01/feb/2016
 *      Author: lorenzo -> dan
 */

#include "PatchyParticleDan.h"
#include "../Utilities/oxDNAException.h"

template<typename number>
PatchyParticleDan<number>::PatchyParticleDan(int _N_patches, LR_vector<number> * _inp_patch_vectors, LR_vector<number> * _inp_ref_vectors, bool tor_flag) : BaseParticle<number>() {
//Original: PatchyParticle<number>::PatchyParticle(char patch_file[], int particle_number) : BaseParticle<number>() {

        //printf("PP, PatchyParticleDan\n");

        //printf("PP 1\n");

        /*for (int i = 0; i < N_patches; i++) {
          printf("patch %d patch vector %f %f %f\n", i, patch_vector[i].x, patch_vector[i].y, patch_vector[i].z);
        }*/

        //tor_flag = flag;
	/*if (_N_patches == 4) {
	  printf("PP tor_flag %d\n", tor_flag);
	  }*/

        _base_patch_vectors = NULL;
        _base_ref_vectors = NULL;
        _ref_vectors = NULL;

	this->N_int_centers = _N_patches;
	//printf("PP this->N_int_centers %d\n", this->N_int_centers);
        //printf("PP 2\n");

	//Create new arrays (one row for each patch)
	this->int_centers = new LR_vector<number>[_N_patches];
        //printf("PP 3\n");
	_base_patch_vectors = new LR_vector<number>[_N_patches];
        //printf("PP 4\n");
	_base_ref_vectors = new LR_vector<number>[_N_patches];
	_ref_vectors = new LR_vector<number>[_N_patches];
	//printf("PP 5\n");

	//Define _base_patch_vectors and _base_ref_vectors
	for(int patch = 0; patch < this->N_int_centers; patch++) {

	  _base_patch_vectors[patch] = _inp_patch_vectors[patch];
	  _base_ref_vectors[patch] = _inp_ref_vectors[patch];
	  //printf("PP 6 patch %d _base_patch_vectors %f %f %f\n", patch, _base_patch_vectors[patch].x, _base_patch_vectors[patch].y, _base_patch_vectors[patch].z);

	}

	//printf("PP 7\n");
	
	/*for(int patch = 0; patch < this->N_int_centers; patch++) {
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

template<typename number>
PatchyParticleDan<number>::~PatchyParticleDan() {
        //printf("PP, ~PatchyParticleDan\n");

        /*Not sure if this is needed, and if it is, should it follow format below? 10/11/16
	  delete[] this->int_centers;*/
	if(_base_patch_vectors != NULL) delete[] _base_patch_vectors;
	if(_base_ref_vectors != NULL) delete[] _base_ref_vectors;
	if(_ref_vectors != NULL) delete[] _ref_vectors;

}

template<typename number>
void PatchyParticleDan<number>::set_positions() {
        //printf("PP, set_positions\n");

        //Set patch vectors and reference vectors based on particle orientation
        for(int patch = 0; patch < this->N_int_centers; patch++) {

	  this->int_centers[patch] = this->orientation * _base_patch_vectors[patch];
	  _ref_vectors[patch] = this->orientation * _base_ref_vectors[patch];

	}

        //printf("PP 10\n");

}

template<typename number>
void PatchyParticleDan<number>::copy_from(const PatchyParticleDan<number> &p) {

        //To inherit/extend copy_from from BaseParticle (see Lorenzo email) 
        BaseParticle<number>::copy_from(p);
	//To add the same functionality for ref_vectors
	for(int patch = 0; patch < this->N_int_centers; patch++) _ref_vectors[patch] = p._ref_vectors[patch];

}

template class PatchyParticleDan<float>;
template class PatchyParticleDan<double>;
