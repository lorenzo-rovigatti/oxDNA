/*
 * SPBAnalysis.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#include "SPBAnalysis.h"
#include <sstream>

template<typename number>
SPBAnalysis<number>::SPBAnalysis(): BaseObservable<number>() {
	_N_bins = -1;
	_confs = 0;
}

template<typename number>
SPBAnalysis<number>::~SPBAnalysis() {

}


template<typename number>
void SPBAnalysis<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable<number>::get_settings(my_inp, sim_inp);

	getInputNumber(&my_inp, "bin", &_bin, 1);
}

template<typename number>
void SPBAnalysis<number>::init(ConfigInfo<number> &config_info) {
	BaseObservable<number>::init(config_info);

	_N_bins = ceil(*config_info.box_side / _bin / 2.);
	_cx.resize(_N_bins, 0.);
	_cy.resize(_N_bins, 0.);
	_cz.resize(_N_bins, 0.);
}

template<typename number>
std::string SPBAnalysis<number>::get_output_string(llint curr_step) {
	LR_vector<number> c_pos = this->_config_info.particles[0]->pos;

	for(int i = 1; i < *this->_config_info.N; i++) {
		LR_vector<number> rel_pos = c_pos - this->_config_info.particles[i]->pos;
		_cx[floor(fabs(rel_pos.x) / _bin)]++;
		_cy[floor(fabs(rel_pos.y) / _bin)]++;
		_cz[floor(fabs(rel_pos.z) / _bin)]++;
	}
	_confs++;

	stringstream ss;
	for(int i = 0; i < _N_bins; i++) {
		number factor = _bin * _confs;
		ss << (i + 0.5) * _bin << " ";
		ss << _cx[i] / factor << " ";
		ss << _cy[i] / factor << " ";
		ss << _cz[i] / factor << endl;
	}

	return ss.str();
}

template class SPBAnalysis<float>;
template class SPBAnalysis<double>;
