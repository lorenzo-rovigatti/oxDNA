/*
 * CubicBox.cpp
 *
 *  Created on: 17/mar/2015
 *      Author: lorenzo
 */

#include "CubicBox.h"

#include "../Utilities/Utils.h"
#include "../Particles/BaseParticle.h"

#include <limits>


using namespace std;

CubicBox::CubicBox() {
	_side = -1.0;
}

CubicBox::~CubicBox() {

}

void CubicBox::get_settings(input_file &inp) {

}

void CubicBox::init(number Lx, number Ly, number Lz) {
	if(Lx != Ly || Ly != Lz || Lz != Lx) throw oxDNAException("The box in the configuration file is not cubic (%f %f %f). Non-cubic boxes can be used by adding a 'box_type = orthogonal' option.", Lx, Ly, Lz);

	_side = Lx;
	_sides.x = _sides.y = _sides.z = Lx;

	CONFIG_INFO->notify(INIT_EVENT);
}

LR_vector CubicBox::normalised_in_box(const LR_vector &v) {
	return LR_vector(
		(v.x / _side - floor(v.x / _side)) * (1.f - std::numeric_limits<number>::epsilon()),
		(v.y / _side - floor(v.y / _side)) * (1.f - std::numeric_limits<number>::epsilon()),
		(v.z / _side - floor(v.z / _side)) * (1.f - std::numeric_limits<number>::epsilon())
	);
}

LR_vector CubicBox::box_sides() const {
	return _sides;
}

inline LR_vector CubicBox::min_image(const LR_vector &v1, const LR_vector &v2) const {
	return LR_vector(
		v2.x - v1.x - rint((v2.x - v1.x) / _side) * _side,
		v2.y - v1.y - rint((v2.y - v1.y) / _side) * _side,
		v2.z - v1.z - rint((v2.z - v1.z) / _side) * _side
	);
}

number CubicBox::sqr_min_image_distance(const LR_vector &v1, const LR_vector &v2) const {
	return min_image(v1, v2).norm();
}

LR_vector CubicBox::get_abs_pos(BaseParticle *p) {
	return p->pos + LR_vector (
			_side * (number)p->_pos_shift[0],
			_side * (number)p->_pos_shift[1],
			_side * (number)p->_pos_shift[2]);
}

void CubicBox::shift_particle(BaseParticle * p, LR_vector &amount) {
	p->_pos_shift[0] += (int) floor(amount.x / _side);
	p->_pos_shift[1] += (int) floor(amount.y / _side);
	p->_pos_shift[2] += (int) floor(amount.z / _side);
	p->pos.x -= _side * floor(amount.x / _side);
	p->pos.y -= _side * floor(amount.y / _side);
	p->pos.z -= _side * floor(amount.z / _side);
}
