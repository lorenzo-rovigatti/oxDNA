/*
 * COMForce.cpp
 *
 *  Created on: 16 May 2014
 *      Author: lorenzo
 */

#include "COMForce.h"

#include "../Utilities/Utils.h"
#include "../Boxes/BaseBox.h"

#include <vector>

using namespace std;

COMForce::COMForce() {
	_rate = 0;
}

COMForce::~COMForce() {

}

std::tuple<std::vector<int>, std::string> COMForce::init(input_file &inp) {
	BaseForce::init(inp);

	getInputString(&inp, "com_list", _com_string, 1);
	getInputString(&inp, "ref_list", _ref_string, 1);
	getInputNumber(&inp, "stiff", &_stiff, 1);
	getInputNumber(&inp, "r0", &_r0, 1);
	getInputNumber(&inp, "rate", &_rate, 0);

	auto com_indexes = Utils::get_particles_from_string(CONFIG_INFO->particles(), _com_string, "COMForce");
	for(auto it = com_indexes.begin(); it != com_indexes.end(); it++) {
		_com_list.insert(CONFIG_INFO->particles()[*it]);
	}

	auto ref_indexes = Utils::get_particles_from_string(CONFIG_INFO->particles(), _ref_string, "COMForce");
	for(auto it = ref_indexes.begin(); it != ref_indexes.end(); it++) {
		_ref_list.insert(CONFIG_INFO->particles()[*it]);
	}

	std::string description = Utils::sformat("COM force of stiffness = %lf, r0 = %lf, rate = %f", _stiff, _r0, _rate);
	return std::make_tuple(com_indexes, description);
}

void COMForce::_compute_coms(llint step) {
	if(step != _last_step) {
		_com = _ref_com = LR_vector(0, 0, 0);
		for(auto p : _com_list) {
			_com += CONFIG_INFO->box->get_abs_pos(p);
		}
		_com /= _com_list.size();

		for(auto p : _ref_list) {
			_ref_com += CONFIG_INFO->box->get_abs_pos(p);
		}
		_ref_com /= _ref_list.size();

		_last_step = step;
	}
}

LR_vector COMForce::value(llint step, LR_vector &pos) {
	_compute_coms(step);
	LR_vector dist = (_ref_com - _com);
	number d_com = dist.module();
	number force = (d_com - (_r0 + _rate * step)) * _stiff / _com_list.size();

	return dist * (force / d_com);
}

number COMForce::potential(llint step, LR_vector &pos) {
	_compute_coms(step);
	return 0.5 * _stiff * SQR((_ref_com - _com).module() - (_r0 + _rate * step)) / _com_list.size();
}
