/*
 * SPBAnalysis.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#include "SPBAnalysis.h"
#include <sstream>

SPBAnalysis::SPBAnalysis() :
				BaseObservable() {
	_N_bins = -1;
	_confs = 0;
}

SPBAnalysis::~SPBAnalysis() {

}

void SPBAnalysis::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputNumber(&my_inp, "bin", &_bin, 1);

	CHECK_BOX("SPBAnalysis", my_inp);
}

void SPBAnalysis::init() {
	BaseObservable::init();

	number box_side = _config_info->box->box_sides().x;
	_N_bins = ceil(box_side / _bin / 2.);
	_cx.resize(_N_bins, 0.);
	_cy.resize(_N_bins, 0.);
	_cz.resize(_N_bins, 0.);
}

std::string SPBAnalysis::get_output_string(llint curr_step) {
	LR_vector c_pos = _config_info->particles()[0]->pos;

	for(int i = 1; i < _config_info->N(); i++) {
		LR_vector rel_pos = c_pos - _config_info->particles()[i]->pos;
		_cx[floor(fabs(rel_pos.x) / _bin)]++;
		_cy[floor(fabs(rel_pos.y) / _bin)]++;
		_cz[floor(fabs(rel_pos.z) / _bin)]++;
	}
	_confs++;

	std::stringstream ss;
	for(int i = 0; i < _N_bins; i++) {
		number factor = _bin * _confs;
		ss << (i + 0.5) * _bin << " ";
		ss << _cx[i] / factor << " ";
		ss << _cy[i] / factor << " ";
		ss << _cz[i] / factor << std::endl;
	}

	return ss.str();
}
