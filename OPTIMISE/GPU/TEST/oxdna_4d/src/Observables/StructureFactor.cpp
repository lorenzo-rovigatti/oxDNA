/*
 * StructureFactor.cpp
 *
 *  Created on: Feb 15, 2016
 *      Author: Lorenzo
 */

#include "StructureFactor.h"

#include <algorithm>

StructureFactor::StructureFactor() {
	_type = -1;
	_max_qs_in_interval = 30;
	_max_qs_delta = 0.001;
	_always_reset = false;
	_max_q = 10;
}

StructureFactor::~StructureFactor() {

}

void StructureFactor::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputNumber(&my_inp, "max_q", &_max_q, 1);
	getInputInt(&my_inp, "int_type", &_type, 0);
	getInputInt(&my_inp, "max_qs_in_interval", &_max_qs_in_interval, 0);
	getInputNumber(&my_inp, "max_qs_delta", &_max_qs_delta, 0);
	getInputBool(&my_inp, "always_reset", &_always_reset, 0);
}

struct sort_qs {
	inline bool operator()(const LR_vector& v1, const LR_vector& v2) {
		return (v1.norm() < v2.norm());
	}
};

void StructureFactor::init() {
	BaseObservable::init();

	LR_vector box_sides = _config_info->box->box_sides();
	number sqr_max_q = SQR(_max_q);
	LR_vector delta_q(2. * M_PI / box_sides.x, 2. * M_PI / box_sides.y, 2. * M_PI / box_sides.z);
	for(int nx = 0; nx <= _max_q / delta_q.x; nx++) {
		for(int ny = 0; ny <= _max_q / delta_q.y; ny++) {
			for(int nz = 0; nz <= _max_q / delta_q.z; nz++) {
				if(nx == 0 && ny == 0 && nz == 0) continue;

				LR_vector new_q(delta_q);
				new_q.x *= nx;
				new_q.y *= ny;
				new_q.z *= nz;

				if(new_q.norm() <= sqr_max_q) _qs.push_back(new_q);
			}
		}
	}

	_qs.sort(sort_qs());

	int q_count = 0;
	double first_q = -1.;
	for(typename std::list<LR_vector >::iterator it = _qs.begin(); it != _qs.end();) {
		q_count++;
		double q_mod = it->module();
		if(first_q < 0.) first_q = q_mod;

		if(q_count > _max_qs_in_interval) it = _qs.erase(it);
		else it++;

		if(it == _qs.end() || fabs(it->norm() - SQR(first_q)) > _max_qs_delta) {
			q_count = 0;
			first_q = -1.;
		}
	}

	_sq.resize(_qs.size(), 0.);
	_sq_cos.resize(_qs.size(), 0.);
	_sq_sin.resize(_qs.size(), 0.);

	OX_LOG(Logger::LOG_INFO, "StructureFactor: %d wave vectors", _qs.size());
}

void StructureFactor::update_data(llint curr_step) {
	int N = _config_info->N();
	uint32_t nq = 0;
	for(typename std::list<LR_vector >::iterator it = _qs.begin(); it != _qs.end(); nq++, it++) {
		double sq_cos = 0.;
		double sq_sin = 0.;
		int N_type = 0;
		for(int i = 0; i < N; i++) {
			BaseParticle *p = _config_info->particles()[i];
			if(_type == -1 || p->type == _type) {
				LR_vector r(p->pos.x, p->pos.y, p->pos.z);
				number qr = *it * r;
				sq_cos += cos(qr);
				sq_sin += sin(qr);
				N_type++;
			}
		}

		_sq[nq] += (SQR(sq_cos) + SQR(sq_sin)) / N_type;
	}
	_times_updated++;
}

std::string StructureFactor::get_output_string(llint curr_step) {
	if(_update_every == 0) {
		update_data(curr_step);
	}

	std::stringstream ret;
	ret.precision(9);
	int q_count = 0;
	double avg_q_mod = 0.;
	double sq_mean = 0.;
	double first_q = -1;
	uint32_t nq = 0;
	for(auto it = _qs.begin(); it != _qs.end(); nq++) {
		q_count++;
		double q_mod = it->module();
		if(first_q < 0.) first_q = q_mod;
		avg_q_mod += q_mod;
		sq_mean += _sq[nq];
		it++;
		if(nq + 1 == _qs.size() || fabs(it->norm() - SQR(first_q)) > _max_qs_delta) {
			avg_q_mod /= q_count;
			sq_mean /= q_count * _times_updated;
			ret << avg_q_mod << " " << sq_mean << std::endl;
			q_count = 0;
			avg_q_mod = sq_mean = 0.;
			first_q = -1.;
		}
	}

	if(_always_reset) {
		_times_updated = 0;
		std::fill(_sq.begin(), _sq.end(), 0.);
	}

	return ret.str();
}
