/*
 * CubicBox.cpp
 *
 *  Created on: 17/mar/2015
 *      Author: lorenzo
 */

#include "CubicBox.h"

#include <cfloat>
#include "../Utilities/Utils.h"

using namespace std;

template<typename number>
CubicBox<number>::CubicBox() {

}

template<typename number>
CubicBox<number>::~CubicBox() {

}

template<typename number>
void CubicBox<number>::get_settings(input_file &inp) {

}

template<typename number>
void CubicBox<number>::init(number Lx, number Ly, number Lz) {
	if(Lx != Ly || Ly != Lz || Lz != Lx) throw oxDNAException("The box in the configuration file is not cubic (%f %f %f)", Lx, Ly, Lz);

	_side = Lx;
	_sides.x = _sides.y = _sides.z = Lx;
}

template<>
LR_vector<float> CubicBox<float>::normalised_in_box(const LR_vector<float> &v) {
	return LR_vector<float> (
		(v.x / _side - floorf(v.x / _side)) * (1.f - FLT_EPSILON),
		(v.y / _side - floorf(v.y / _side)) * (1.f - FLT_EPSILON),
		(v.z / _side - floorf(v.z / _side)) * (1.f - FLT_EPSILON)
	);
}

template<>
LR_vector<double> CubicBox<double>::normalised_in_box(const LR_vector<double> &v) {
	return LR_vector<double> (
		(v.x / _side - floorf(v.x / _side)) * (1.f - DBL_EPSILON),
		(v.y / _side - floorf(v.y / _side)) * (1.f - DBL_EPSILON),
		(v.z / _side - floorf(v.z / _side)) * (1.f - DBL_EPSILON)
	);
}

template<typename number>
LR_vector<number> &CubicBox<number>::box_sides()  {
	return _sides;
}

template<typename number>
LR_vector<number> CubicBox<number>::min_image(const LR_vector<number> &v1, const LR_vector<number> &v2) const {
	return LR_vector<number> (
		v2.x - v1.x - rint((v2.x - v1.x) / _side) * _side,
		v2.y - v1.y - rint((v2.y - v1.y) / _side) * _side,
		v2.z - v1.z - rint((v2.z - v1.z) / _side) * _side
	);
}

template<typename number>
number CubicBox<number>::sqr_min_image_distance(const LR_vector<number> &v1, const LR_vector<number> &v2) const {
	number nx = v2.x - v1.x;
	number ny = v2.y - v1.y;
	number nz = v2.z - v1.z;

	nx -= rint(nx / _side) * _side;
	ny -= rint(ny / _side) * _side;
	nz -= rint(nz / _side) * _side;

	return nx*nx + ny*ny + nz*nz;
}

template<typename number>
void CubicBox<number>::apply_boundary_conditions(BaseParticle<number> **particles, int N) {

}

template class CubicBox<float>;
template class CubicBox<double>;
