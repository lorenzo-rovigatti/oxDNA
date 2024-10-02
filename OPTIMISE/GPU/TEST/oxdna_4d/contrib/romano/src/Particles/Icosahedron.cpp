/*
 * Icosahedron.cpp
 *
 *  Created on: 25/lug/2018
 *      Author: Flavio
 */

#include "Icosahedron.h"
#include "Utilities/oxDNAException.h"

template<typename number>
Icosahedron<number>::Icosahedron(int N_patches) : BaseParticle<number>() {
	this->N_int_centers = 12;
	this->int_centers = new LR_vector<number>[12];
	_N_patches = N_patches;
	
	_vertexes = new LR_vector<number>[12];
	_patch_indexes = new int[N_patches];

	_set_vertexes();
}

template<typename number>
Icosahedron<number>::~Icosahedron() {
	delete[] _vertexes;
	delete[] _patch_indexes;
}

template<typename number>
void Icosahedron<number>::_set_vertexes() {
	double t = (1. + sqrt(5.))/2;      // golden radius
	double r = 2. * sin(2. * M_PI/5);  // radius of icosahedron of edge lenght two
	number a = 1./ (2. * r);           // will generate patches at a distance 0.5 from center
	number b = t / (2. * r);

	// 12 vertexes of the icosahedron; circular
	// permutations of (0, +/-a, +/-b)
	// each vertex +6 is the opposite; this is 
	// helpful in the interaction
	_vertexes[ 0] = LR_vector<number>( 0.,  a,  b); //   0
	_vertexes[ 1] = LR_vector<number>( 0.,  a, -b); //   1
	_vertexes[ 2] = LR_vector<number>(  b, 0.,  a); //   4
	_vertexes[ 3] = LR_vector<number>(  b, 0., -a); //   5
	_vertexes[ 4] = LR_vector<number>(  a,  b, 0.); //   8
	_vertexes[ 5] = LR_vector<number>(  a, -b, 0.); //   9
	_vertexes[ 6] = LR_vector<number>( 0., -a, -b); //  -0
	_vertexes[ 7] = LR_vector<number>( 0., -a,  b); //  -1
	_vertexes[ 8] = LR_vector<number>( -b, 0., -a); //  -4
	_vertexes[ 9] = LR_vector<number>( -b, 0.,  a); //  -5
	_vertexes[10] = LR_vector<number>( -a, -b, 0.); //  -8
	_vertexes[11] = LR_vector<number>( -a,  b, 0.); //  -9

	switch(this->_N_patches) {
	case 1: {
		_patch_indexes[0] = 0;
		break;
	}
	case 2: {
		_patch_indexes[0] = 0;
		_patch_indexes[1] = 3;
		break;
	}
	case 6: {
		_patch_indexes[0] = 0;
		_patch_indexes[1] = 2;
		_patch_indexes[2] = 6;
		_patch_indexes[3] = 1;
		_patch_indexes[4] = 3;
		_patch_indexes[5] = 5;
		break;
	}
	default:
		throw oxDNAException("Unsupported number of patches %d\n", this->N_int_centers);
	}

	this->set_positions();
}

template<typename number>
void Icosahedron<number>::set_positions() {
	for(int i = 0; i < this->N_int_centers; i++) this->int_centers[i] = (this->orientation*_vertexes[i]);
}

template class Icosahedron<float>;
template class Icosahedron<double>;
