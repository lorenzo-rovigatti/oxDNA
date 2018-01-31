/*
 * LeesEdwardsCubicBox.cpp
 *
 *  Created on: 18 Aug 2017
 *      Author: lorenzo
 */

#include "LeesEdwardsCubicBox.h"

#include "../Utilities/ConfigInfo.h"

template<typename number>
LeesEdwardsCubicBox<number>::LeesEdwardsCubicBox() {

}

template<typename number>
LeesEdwardsCubicBox<number>::~LeesEdwardsCubicBox() {

}

template<typename number>
void LeesEdwardsCubicBox<number>::get_settings(input_file &inp) {
	number dt, shear_rate;
	getInputNumber(&inp, "dt", &dt, 1);
	getInputNumber(&inp, "lees_edwards_shear_rate", &shear_rate, 1);
	_factor = dt * shear_rate;
}

template<typename number>
void LeesEdwardsCubicBox<number>::init(number Lx, number Ly, number Lz) {
	CubicBox<number>::init(Lx, Ly, Lz);

	_factor *= Ly;
}

template<typename number>
LR_vector<number> LeesEdwardsCubicBox<number>::min_image(const LR_vector<number> &v1, const LR_vector<number> &v2) const {
	static llint last_step = -1;
	static number delta_x = 0.;

	llint curr_step = CONFIG_INFO->curr_step;
	if(curr_step != last_step) {
		delta_x = _factor*curr_step;;
		last_step = curr_step;
	}

	number ny = v2.y - v1.y;
	number cy = rint(ny / this->_side);
	number nx = v2.x - v1.x - cy * delta_x;

	return LR_vector<number>(
			nx - rint(nx / this->_side) * this->_side,
			ny - cy * this->_side,
			v2.z - v1.z - rint((v2.z - v1.z) / this->_side) * this->_side
	);
}

template<typename number>
number LeesEdwardsCubicBox<number>::sqr_min_image_distance(const LR_vector<number> &v1, const LR_vector<number> &v2) const {
	return min_image(v1, v2).norm();
}

//template<typename number>
//void LeesEdwardsCubicBox<number>::shift_particle (BaseParticle<number> *p, LR_vector<number> &amount) {
//	p->_pos_shift[0] += (int) floor(amount.x / this->_side);
//	p->_pos_shift[2] += (int) floor(amount.z / this->_side);
//	p->pos.x -= this->_side * floor(amount.x / this->_side);
//	p->pos.z -= this->_side * floor(amount.z / this->_side);
//}

template class LeesEdwardsCubicBox<float> ;
template class LeesEdwardsCubicBox<double> ;
