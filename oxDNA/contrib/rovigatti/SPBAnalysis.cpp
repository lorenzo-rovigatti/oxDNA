/*
 * SPBAnalysis.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#include "SPBAnalysis.h"

template<typename number>
SPBAnalysis<number>::SPBAnalysis(): BaseObservable<number>() {
	_N_bins = -1;
	_shift = -1;
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

	_N_bins = ceil(*config_info.box_side / _bin);
	_cx.resize(_N_bins);
	_cy.resize(_N_bins);
	_cz.resize(_N_bins);
	_shift = _N_bins / 2;
}

template<typename number>
std::string SPBAnalysis<number>::get_output_string(llint curr_step) {
	LR_vector<number> c_pos = this->_config_info.particles[0]->pos;

	for(int i = 0; i < _N_bins; i++) _cx[i] = _cy[i] = _cz[i] = 0.;

	for(int i = 1; i < *this->_config_info.N; i++) {
		LR_vector<number> rel_pos = c_pos - this->_config_info.particles[i]->pos;
		_cx[floor(rel_pos.x / _bin) + _shift]++;
		_cy[floor(rel_pos.y / _bin) + _shift]++;
		_cz[floor(rel_pos.z / _bin) + _shift]++;
	}

	number lx = 0.;
	number ly = 0.;
	number lz = 0.;

	for(int i = 0; i < _N_bins; i++) {
		number dist = ((i - _shift) + 0.5) * _bin;
		lx += _cx[i] * SQR(dist) * _bin;
		ly += _cy[i] * SQR(dist) * _bin;
		lz += _cz[i] * SQR(dist) * _bin;
	}

	lx = sqrt(lx);
	ly = sqrt(ly);
	lz = sqrt(lz);

	return Utils::sformat("% 10.6lf % 10.6lf % 10.6lf", lx, ly, lz);
}

template class SPBAnalysis<float>;
template class SPBAnalysis<double>;
