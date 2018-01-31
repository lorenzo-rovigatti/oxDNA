/*
 * FormFactor.cpp
 *
 *  Created on: Apr 21 2017
 *      Author: Lorenzo
 */

#include "FormFactor.h"

#include <algorithm>

template<typename number>
FormFactor<number>::FormFactor() {
	_min_q = _max_q = -1;
	_mult_q = 1.2;
	_type = -1;
	_n_qs = 30;
}

template<typename number>
FormFactor<number>::~FormFactor() {

}

template<typename number>
void FormFactor<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	getInputNumber(&my_inp, "min_q", &_min_q, 1);
	getInputNumber(&my_inp, "max_q", &_max_q, 1);
	getInputNumber(&my_inp, "mult_q", &_mult_q, 0);
	getInputInt(&my_inp, "int_type", &_type, 0);
	getInputInt(&my_inp, "n_qs", &_n_qs, 0);
}

struct sort_qs {
	inline bool operator()(const LR_vector<double>& v1, const LR_vector<double>& v2) {
		return (v1.norm() < v2.norm());
	}
};

template<typename number>
void FormFactor<number>::init(ConfigInfo<number> &config_info) {
	BaseObservable<number>::init(config_info);

	double curr_mod = _min_q;
	int tot_n_qs = 0;
	while(curr_mod < _max_q) {
		curr_mod *= _mult_q;
		tot_n_qs++;
	}
	tot_n_qs *= _n_qs;
	OX_LOG(Logger::LOG_INFO, "FormFactor: %d wave vectors", tot_n_qs);
}

template<typename number>
std::string FormFactor<number>::get_output_string(llint curr_step) {
	stringstream ret;
	ret.precision(9);

	int N = *this->_config_info.N;

	LR_vector<number> com(0., 0., 0.);
	for(int j = 0; j < N; j++) {
		BaseParticle<number> *p = this->_config_info.particles[j];
		com += this->_config_info.box->get_abs_pos(p);
	}
	com /= N;

	double curr_mod = _min_q;
	while(curr_mod < _max_q) {
		double sq_avg = 0.;
		for(int i = 0; i < _n_qs; i++) {
			LR_vector<double> new_q = Utils::get_random_vector<double>()*curr_mod;
			int N_type = 0;
			double sq_cos = 0.;
			double sq_sin = 0.;
			for(int j = 0; j < N; j++) {
				BaseParticle<number> *p = this->_config_info.particles[j];
				if(_type == -1 || p->type == _type) {
					LR_vector<number> my_pos = this->_config_info.box->get_abs_pos(p) - com;
					LR_vector<double> r(my_pos.x, my_pos.y, my_pos.z);
					number qr = new_q*r;
					sq_cos += cos(qr);
					sq_sin += sin(qr);
					N_type++;
				}
			}
			sq_avg += (SQR(sq_cos) + SQR(sq_sin))/N_type;
		}
		sq_avg /= _n_qs;
		ret << curr_mod << " " << sq_avg << endl;

		curr_mod *= _mult_q;
	}

	return ret.str();
}

template class FormFactor<float>;
template class FormFactor<double>;
