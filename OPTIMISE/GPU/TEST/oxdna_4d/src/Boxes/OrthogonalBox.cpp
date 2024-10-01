/*
 * OrthogonalBox.cpp
 *
 *  Created on: 17/mar/2015
 *      Author: lorenzo
 */

#include "OrthogonalBox.h"

#include "../Utilities/Utils.h"
#include "../Particles/BaseParticle.h"

#include <limits>

using namespace std;

OrthogonalBox::OrthogonalBox() :
				BaseBox() {

}

OrthogonalBox::~OrthogonalBox() {

}

void OrthogonalBox::get_settings(input_file &inp) {

}

void OrthogonalBox::init(number Lx, number Ly, number Lz) {
	_sides.x = Lx;
	_sides.y = Ly;
	_sides.z = Lz;

	CONFIG_INFO->notify(INIT_EVENT);
}

LR_vector OrthogonalBox::normalised_in_box(const LR_vector &v) {
	return LR_vector((v.x / _sides.x - floor(v.x / _sides.x)) * (1.f - std::numeric_limits<number>::epsilon()), (v.y / _sides.y - floor(v.y / _sides.y)) * (1.f - std::numeric_limits<number>::epsilon()), (v.z / _sides.z - floor(v.z / _sides.z)) * (1.f - std::numeric_limits<number>::epsilon()));
}

LR_vector OrthogonalBox::box_sides() const {
	return _sides;
}

LR_vector OrthogonalBox::min_image(const LR_vector &v1, const LR_vector &v2) const {
	return LR_vector(v2.x - v1.x - rint((v2.x - v1.x) / _sides.x) * _sides.x, v2.y - v1.y - rint((v2.y - v1.y) / _sides.y) * _sides.y, v2.z - v1.z - rint((v2.z - v1.z) / _sides.z) * _sides.z);
}

number OrthogonalBox::sqr_min_image_distance(const LR_vector &v1, const LR_vector &v2) const {
	number nx = v2.x - v1.x;
	number ny = v2.y - v1.y;
	number nz = v2.z - v1.z;

	nx -= rint(nx / _sides.x) * _sides.x;
	ny -= rint(ny / _sides.y) * _sides.y;
	nz -= rint(nz / _sides.z) * _sides.z;

	return nx * nx + ny * ny + nz * nz;
}

LR_vector OrthogonalBox::get_abs_pos(BaseParticle * p) {
	return p->pos + LR_vector(_sides.x * (number) p->_pos_shift[0], _sides.y * (number) p->_pos_shift[1], _sides.z * (number) p->_pos_shift[2]);
}

void OrthogonalBox::shift_particle(BaseParticle * p, LR_vector &amount) {
	p->_pos_shift[0] += (int) floor(amount.x / _sides.x);
	p->_pos_shift[1] += (int) floor(amount.y / _sides.y);
	p->_pos_shift[2] += (int) floor(amount.z / _sides.z);
	p->pos.x -= _sides.x * floor(amount.x / _sides.x);
	p->pos.y -= _sides.y * floor(amount.y / _sides.y);
	p->pos.z -= _sides.z * floor(amount.z / _sides.z);
}
