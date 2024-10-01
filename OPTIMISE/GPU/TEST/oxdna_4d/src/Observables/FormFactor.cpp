/*
 * FormFactor.cpp
 *
 *  Created on: Apr 21 2017
 *      Author: Lorenzo
 */

#include "FormFactor.h"

#include <algorithm>

FormFactor::FormFactor() {
	_min_q = _max_q = -1;
	_mult_q = 1.2;
	_type = -1;
	_n_qs = 30;
}

FormFactor::~FormFactor() {

}

void FormFactor::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputNumber(&my_inp, "min_q", &_min_q, 1);
	getInputNumber(&my_inp, "max_q", &_max_q, 1);
	getInputNumber(&my_inp, "mult_q", &_mult_q, 0);
	getInputInt(&my_inp, "int_type", &_type, 0);
	getInputInt(&my_inp, "n_qs", &_n_qs, 0);
}

struct sort_qs {
	inline bool operator()(const LR_vector &v1, const LR_vector& v2) {
		return (v1.norm() < v2.norm());
	}
};

void FormFactor::init() {
	BaseObservable::init();

	double curr_mod = _min_q;
	int tot_n_qs = 0;
	while(curr_mod < _max_q) {
		curr_mod *= _mult_q;
		tot_n_qs++;
	}
	tot_n_qs *= _n_qs;
	OX_LOG(Logger::LOG_INFO, "FormFactor: %d wave vectors", tot_n_qs);
}

std::string FormFactor::get_output_string(llint curr_step) {
	std::stringstream ret;
	ret.precision(9);

	int N = _config_info->N();

	LR_vector com(0., 0., 0.);
	for(int j = 0; j < N; j++) {
		BaseParticle *p = _config_info->particles()[j];
		com += _config_info->box->get_abs_pos(p);
	}
	com /= N;

	double curr_mod = _min_q;
	while(curr_mod < _max_q) {
		double sq_avg = 0.;
		for(int i = 0; i < _n_qs; i++) {
			LR_vector new_q = Utils::get_random_vector() * curr_mod;
			int N_type = 0;
			double sq_cos = 0.;
			double sq_sin = 0.;
			for(int j = 0; j < N; j++) {
				BaseParticle *p = _config_info->particles()[j];
				if(_type == -1 || p->type == _type) {
					LR_vector my_pos = _config_info->box->get_abs_pos(p) - com;
					LR_vector r(my_pos.x, my_pos.y, my_pos.z);
					number qr = new_q * r;
					sq_cos += cos(qr);
					sq_sin += sin(qr);
					N_type++;
				}
			}
			sq_avg += (SQR(sq_cos) + SQR(sq_sin)) / N_type;
		}
		sq_avg /= _n_qs;
		ret << curr_mod << " " << sq_avg << std::endl;

		curr_mod *= _mult_q;
	}

	return ret.str();
}
