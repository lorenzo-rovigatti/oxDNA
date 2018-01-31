/*
 * StructureFactor.cpp
 *
 *  Created on: Feb 15, 2016
 *      Author: Lorenzo
 */

#include "StructureFactor.h"

#include <algorithm>

template<typename number>
StructureFactor<number>::StructureFactor() {
	_nconf = 0;
	_type = -1;
	_max_qs_in_interval = 30;
	_max_qs_delta = 0.001;
	_always_reset = false;
}

template<typename number>
StructureFactor<number>::~StructureFactor() {

}

template<typename number>
void StructureFactor<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	getInputNumber(&my_inp, "max_q", &_max_q, 1);
	getInputInt(&my_inp, "int_type", &_type, 0);
	getInputInt(&my_inp, "max_qs_in_interval", &_max_qs_in_interval, 0);
	getInputNumber(&my_inp, "max_qs_delta", &_max_qs_delta, 0);
	getInputBool(&my_inp, "always_reset", &_always_reset, 0);
}

struct sort_qs {
	inline bool operator()(const LR_vector<double>& v1, const LR_vector<double>& v2) {
		return (v1.norm() < v2.norm());
	}
};

template<typename number>
void StructureFactor<number>::init(ConfigInfo<number> &config_info) {
	BaseObservable<number>::init(config_info);

	LR_vector<number> box_sides = config_info.box->box_sides();
	number sqr_max_q = SQR(_max_q);
	LR_vector<double> delta_q(2.*M_PI/box_sides.x, 2.*M_PI/box_sides.y, 2.*M_PI/box_sides.z);
	for(int nx = 0; nx <= _max_q/delta_q.x; nx++) {
		for(int ny = 0; ny <= _max_q/delta_q.y; ny++) {
			for(int nz = 0; nz <= _max_q/delta_q.z; nz++) {
				if(nx == 0 && ny == 0 && nz == 0) continue;

				LR_vector<double> new_q(delta_q);
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
	for(typename std::list<LR_vector<double> >::iterator it = _qs.begin(); it != _qs.end();) {
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

template<typename number>
std::string StructureFactor<number>::get_output_string(llint curr_step) {
	if(_always_reset) {
		_nconf = 1;
		std::fill(_sq.begin(), _sq.end(), 0.);
	}
	else _nconf += 1;

	int N = *this->_config_info.N;
	uint nq = 0;
	for(typename std::list<LR_vector<double> >::iterator it = _qs.begin(); it != _qs.end(); nq++, it++) {
		double sq_cos = 0.;
		double sq_sin = 0.;
		int N_type = 0;
		for(int i = 0; i < N; i++) {
			BaseParticle<number> *p = this->_config_info.particles[i];
			if(_type == -1 || p->type == _type) {
				LR_vector<double> r(p->pos.x, p->pos.y, p->pos.z);
				number qr = *it*r;
				sq_cos += cos(qr);
				sq_sin += sin(qr);
				N_type++;
			}
		}

		_sq[nq] += (SQR(sq_cos) + SQR(sq_sin))/N_type;
	}

	stringstream ret;
	ret.precision(9);
	int q_count = 0;
	double avg_q_mod = 0.;
	double sq_mean = 0.;
	double first_q = -1;
	nq = 0;
	for(typename std::list<LR_vector<double> >::iterator it = _qs.begin(); it != _qs.end(); nq++) {
		q_count++;
		double q_mod = it->module();
		if(first_q < 0.) first_q = q_mod;
		avg_q_mod += q_mod;
		sq_mean += _sq[nq];
		it++;
		if(nq+1 == _qs.size() || fabs(it->norm() - SQR(first_q)) > _max_qs_delta) {
			avg_q_mod /= q_count;
			sq_mean /= q_count*_nconf;
			ret << avg_q_mod << " " << sq_mean << endl;
			q_count = 0;
			avg_q_mod = sq_mean = 0.;
			first_q = -1.;
		}
	}

	return ret.str();
}

template class StructureFactor<float>;
template class StructureFactor<double>;
