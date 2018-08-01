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

	// todo: perhaps remove these useless arrays...
	_faces = new int[20 * 3]; // store an array of indexes of vertexes defining a face
	_close_faces = new int[12 * 5];  // store the indexes of the faces close to each vertex 

	_set_vertexes();
}

template<typename number>
Icosahedron<number>::~Icosahedron() {
	delete[] _vertexes;
	delete[] _faces;
	delete[] _patch_indexes;
	delete[] _close_faces;
}

template<typename number>
void Icosahedron<number>::_set_vertexes() {
	double t = (1. + sqrt(5.))/2;      // golden radius
	double r = 2. * sin(2. * M_PI/5);  // radius of icosahedron of edge lenght two
	number a = 1./ (2. * r);           // will generate patches at a distance 0.5 from center
	number b = t / (2. * r);

	// 12 vertexes of the icosahedron; circular
	// permutations of (0, +/-a, +/-b)
	_vertexes[ 0] = LR_vector<number>( 0.,  a,  b);
	_vertexes[ 1] = LR_vector<number>( 0.,  a, -b);
	_vertexes[ 2] = LR_vector<number>( 0., -a,  b);
	_vertexes[ 3] = LR_vector<number>( 0., -a, -b);

	_vertexes[ 4] = LR_vector<number>(  b, 0.,  a);
	_vertexes[ 5] = LR_vector<number>(  b, 0., -a);
	_vertexes[ 6] = LR_vector<number>( -b, 0.,  a);
	_vertexes[ 7] = LR_vector<number>( -b, 0., -a);

	_vertexes[ 8] = LR_vector<number>(  a,  b, 0.);
	_vertexes[ 9] = LR_vector<number>(  a, -b, 0.);
	_vertexes[10] = LR_vector<number>( -a,  b, 0.);
	_vertexes[11] = LR_vector<number>( -a, -b, 0.);

	// we now need to figure out which vertexes belong to which face
	// angle between vertexes of the same face: 1.1071487177940436, ~63.435 degrees
	// the possible angles between all pairs are: ~63.4, ~116.6 e 180.
	// each vertex has 5 vertexes at 63, 5 more at 116 and 1 at 180 (opposite)
	int nface = 0;
	number thres = 0.;  // threshold angle is 90 degrees
	for (int i = 0; i < 12; i ++) {
		for (int j = 0; j < i; j ++) {
			for (int k = 0; k < j; k ++) {
				if ((_vertexes[i]*_vertexes[j] > thres) &&
				    (_vertexes[i]*_vertexes[k] > thres) &&
				    (_vertexes[j]*_vertexes[k] > thres)) {
						_faces[3*nface + 0] = i;
						_faces[3*nface + 1] = j;
						_faces[3*nface + 2] = k;
						/*printf ("\n\n%d %d %d @\n", i, j, k);
						printf ("%7.5g %7.5g %7.5g\n", _vertexes[i].x, _vertexes[i].y, _vertexes[i].z);
						printf ("%7.5g %7.5g %7.5g\n", _vertexes[j].x, _vertexes[j].y, _vertexes[j].z);
						printf ("%7.5g %7.5g %7.5g\n", _vertexes[k].x, _vertexes[k].y, _vertexes[k].z);
						printf ("  %g %g %g\n", 4.f*(_vertexes[i]*_vertexes[j]), 4.f*(_vertexes[i]*_vertexes[k]), 4.*(_vertexes[j]*_vertexes[k]));*/
						nface ++;
				}
			}
		}
	}
	if (nface != 20) throw oxDNAException ("no can do here %d...", nface);

	// we also want an array of faces close to each vertex
	for (int i = 0; i < 12; i ++) {
		int mynfaces = 0;
		for (int j = 0; j < 20; j ++) {
			if ( (_faces[3*j + 0] == i) ||
				 (_faces[3*j + 1] == i) ||
				 (_faces[3*j + 2] == i) ) {
				_close_faces[5 * i + mynfaces] = j;
				mynfaces ++;
			}
		}
		if (mynfaces != 5) throw oxDNAException ("no can do 2");
	}

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
