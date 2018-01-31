/*
 * OrthogonalBox.cpp
 *
 *  Created on: 17/mar/2015
 *      Author: lorenzo
 */

#include "OrthogonalBox.h"

#include <cfloat>
#include "../Utilities/Utils.h"
#include "../Particles/BaseParticle.h"

using namespace std;

template<typename number>
OrthogonalBox<number>::OrthogonalBox() {

}

template<typename number>
OrthogonalBox<number>::~OrthogonalBox() {

}

template<typename number>
void OrthogonalBox<number>::get_settings(input_file &inp) {

}

template<typename number>
void OrthogonalBox<number>::init(number Lx, number Ly, number Lz) {
	_sides.x = Lx;
	_sides.y = Ly;
	_sides.z = Lz;
}

template<>
LR_vector<float> OrthogonalBox<float>::normalised_in_box(const LR_vector<float> &v) {
	return LR_vector<float> (
		(v.x / _sides.x - floorf(v.x / _sides.x)) * (1.f - FLT_EPSILON),
		(v.y / _sides.y - floorf(v.y / _sides.y)) * (1.f - FLT_EPSILON),
		(v.z / _sides.z - floorf(v.z / _sides.z)) * (1.f - FLT_EPSILON)
	);
}

template<>
LR_vector<double> OrthogonalBox<double>::normalised_in_box(const LR_vector<double> &v) {
	return LR_vector<double> (
		(v.x / _sides.x - floorf(v.x / _sides.x)) * (1.f - DBL_EPSILON),
		(v.y / _sides.y - floorf(v.y / _sides.y)) * (1.f - DBL_EPSILON),
		(v.z / _sides.z - floorf(v.z / _sides.z)) * (1.f - DBL_EPSILON)
	);
}

template<typename number>
LR_vector<number> &OrthogonalBox<number>::box_sides()  {
	return _sides;
}

template<typename number>
LR_vector<number> OrthogonalBox<number>::min_image(const LR_vector<number> &v1, const LR_vector<number> &v2) const {
	return LR_vector<number> (
		v2.x - v1.x - rint((v2.x - v1.x) / _sides.x) * _sides.x,
		v2.y - v1.y - rint((v2.y - v1.y) / _sides.y) * _sides.y,
		v2.z - v1.z - rint((v2.z - v1.z) / _sides.z) * _sides.z
	);
}

template<typename number>
number OrthogonalBox<number>::sqr_min_image_distance(const LR_vector<number> &v1, const LR_vector<number> &v2) const {
	number nx = v2.x - v1.x;
	number ny = v2.y - v1.y;
	number nz = v2.z - v1.z;

	nx -= rint(nx / _sides.x) * _sides.x;
	ny -= rint(ny / _sides.y) * _sides.y;
	nz -= rint(nz / _sides.z) * _sides.z;

	return nx*nx + ny*ny + nz*nz;
}

template<typename number>
void OrthogonalBox<number>::apply_boundary_conditions(BaseParticle<number> **particles, int N) {

}

template<typename number>
LR_vector<number> OrthogonalBox<number>::get_abs_pos(BaseParticle<number> * p) {
	return p->pos + LR_vector<number> (
			_sides.x * (number)p->_pos_shift[0],
			_sides.y * (number)p->_pos_shift[1],
			_sides.z * (number)p->_pos_shift[2]);
}

template<typename number>
void OrthogonalBox<number>::shift_particle(BaseParticle<number> * p, LR_vector<number> &amount) {
	p->_pos_shift[0] += (int) floor(amount.x / _sides.x);
	p->_pos_shift[1] += (int) floor(amount.y / _sides.y);
	p->_pos_shift[2] += (int) floor(amount.z / _sides.z);
	p->pos.x -= _sides.x * floor(amount.x / _sides.x);
	p->pos.y -= _sides.y * floor(amount.y / _sides.y);
	p->pos.z -= _sides.z * floor(amount.z / _sides.z);
}

template class OrthogonalBox<float>;
template class OrthogonalBox<double>;
