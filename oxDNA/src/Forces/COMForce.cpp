/*
 * COMForce.cpp
 *
 *  Created on: 16 May 2014
 *      Author: lorenzo
 */

#include <vector>

#include "COMForce.h"

#include "../Utilities/Utils.h"

#include "../Boxes/BaseBox.h"

using namespace std;

template<typename number>
COMForce<number>::COMForce() {
	_r0 = 0;
	_last_step = -1;
	_box_ptr = NULL;
}

template<typename number>
COMForce<number>::~COMForce() {

}

template<typename number>
void COMForce<number>::get_settings(input_file &inp) {
	getInputString(&inp, "com_list", _com_string, 1);
	getInputString(&inp, "ref_list", _ref_string, 1);

	double stiff;
	getInputDouble(&inp, "stiff", &stiff, 1);
	this->_stiff = stiff;

	double r0;
	getInputDouble(&inp, "r0", &r0, 1);
	_r0 = r0;
}

template<typename number>
void COMForce<number>::_check_index(int idx, int N) {
	if(idx < 0 || idx >= N) throw oxDNAException("COMForce: invalid id %d", idx);
}

template<typename number>
void COMForce<number>::init(BaseParticle<number> **particles, int N, BaseBox<number> * box_ptr) {
	_box_ptr = box_ptr;

	vector<int> com_indexes = Utils::getParticlesFromString(particles, N, _com_string, "COMForce");
	for(vector<int>::iterator it = com_indexes.begin(); it != com_indexes.end(); it++) {
		_check_index(*it, N);
		_com_list.insert(particles[*it]);
		particles[*it]->add_ext_force(this);
	}

	vector<int> ref_indexes = Utils::getParticlesFromString(particles, N, _ref_string, "COMForce");
	for(vector<int>::iterator it = ref_indexes.begin(); it != ref_indexes.end(); it++) {
		_check_index(*it, N);
		_ref_list.insert(particles[*it]);
	}

	OX_LOG(Logger::LOG_INFO, "Adding a COM force of stiffness = %lf and r0 = %lf", this->_stiff, _r0);
}

template<typename number>
void COMForce<number>::_compute_coms(llint step) {
	if(step != _last_step) {
		_com = _ref_com = LR_vector<number>(0, 0, 0);
		for(typename set<BaseParticle<number> *>::iterator it = _com_list.begin(); it != _com_list.end(); it++) {
			_com += _box_ptr->get_abs_pos(*it);
		}
		_com /= _com_list.size();

		for(typename set<BaseParticle<number> *>::iterator it = _ref_list.begin(); it != _ref_list.end(); it++) {
			_ref_com += _box_ptr->get_abs_pos(*it);
		}
		_ref_com /= _ref_list.size();

		_last_step = step;
	}
}

template<typename number>
LR_vector<number> COMForce<number>::value(llint step, LR_vector<number> &pos) {
	_compute_coms(step);
	LR_vector<number> dist = (_ref_com - _com);
	number d_com = dist.module();
	number force = (d_com - _r0) * this->_stiff / _com_list.size();

	return dist * (force/d_com);
}

template<typename number>
number COMForce<number>::potential(llint step, LR_vector<number> &pos) {
	_compute_coms(step);
	return 0.5 * this->_stiff * SQR((_ref_com - _com).module() - _r0) / _com_list.size();
}

template class COMForce<double>;
template class COMForce<float>;
