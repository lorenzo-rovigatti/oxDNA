/*
 * GrByInsertion.cpp
 *
 *  Created on: Oct 30, 2013
 *      Author: Lorenzo
 */

#include <sstream>

#include "GrByInsertion.h"

#include "../Interactions/TSPInteraction.h"
#include "Utilities/Utils.h"

using namespace std;

GrByInsertion::GrByInsertion() {
	_n_bins = 0;
	_n_conf = 0;
	_bin = 0.1;
	_T = 0.;
	_inter_hist[AA] = _inter_hist[AB] = _inter_hist[BB] = NULL;
	_inter_norm[AA] = _inter_norm[AB] = _inter_norm[BB] = 0.;
	_insertions = 1000;
	_type = AA;
	_fixed = _movable = P_A;
	_min = _max = 0.;
}

GrByInsertion::~GrByInsertion() {
	if(_inter_hist[AA] != NULL) {
		delete[] _inter_hist[AA];
		delete[] _inter_hist[AB];
		delete[] _inter_hist[BB];
	}
}

void GrByInsertion::get_settings(input_file &my_inp, input_file &sim_inp) {
	float tmp;
	if(getInputFloat(&my_inp, "bin", &tmp, 0) == KEY_FOUND) _bin = tmp;

	getInputInt(&my_inp, "insertions", &_insertions, 0);
	getInputFloat(&sim_inp, "T", &tmp, 1);
	_T = tmp;

	if(getInputFloat(&my_inp, "gr_min", &tmp, 0) == KEY_FOUND) _min = tmp;
	if(getInputFloat(&my_inp, "gr_max", &tmp, 0) == KEY_FOUND) _max = tmp;

	char stype[512];
	getInputString(&my_inp, "gr_type", stype, 1);
	if(!strncmp(stype, "AA", 512)) _type = AA;
	else if(!strncmp(stype, "AB", 512)) _type = AB;
	else if(!strncmp(stype, "BB", 512)) _type = BB;
	else throw oxDNAException("gr_type should be AA, AB, or BB. %s not supported", stype);

	CHECK_BOX("GrByInsertion", my_inp);
}

void GrByInsertion::init() {
	BaseObservable::init();

	if(_max == 0.) _max = _config_info->box->box_sides().x * 0.5;
	_n_bins = (_max - _min) / _bin;
	// rounding
	_bin = (_max - _min) / _n_bins;

	_inter_hist[AA] = new number[_n_bins];
	_inter_hist[AB] = new number[_n_bins];
	_inter_hist[BB] = new number[_n_bins];

	for(int i = 0; i < _n_bins; i++) {
		_inter_hist[AA][i] = 0.;
		_inter_hist[AB][i] = 0.;
		_inter_hist[BB][i] = 0.;
	}

	for(auto p: _config_info->particles()) {
		if(p->strand_id > 1) throw oxDNAException("The system should contain only two chains");
		if(p->type != P_A && p->type != P_B) throw oxDNAException("Only particles of type A and B are allowed");
		_particles[p->strand_id][p->type].push_back(p);
	}

	switch(_type) {
	case AA:
		_fixed = P_A;
		_movable = P_A;
		break;
	case AB:
		_fixed = P_A;
		_movable = P_B;
		break;
	case BB:
		_fixed = P_B;
		_movable = P_B;
		break;
	}
}

int GrByInsertion::_get_bin(number sqr_dist) {
	int bin = (sqrt(sqr_dist) - _min) / _bin;
	if(bin >= _n_bins) bin = -1;

	return bin;
}

LR_vector GrByInsertion::_get_com(int chain, int type) {
	LR_vector res(0., 0., 0.);

	typename vector<BaseParticle *>::iterator it;
	for(it = _particles[chain][type].begin(); it != _particles[chain][type].end(); it++) {
		res = _config_info->box->get_abs_pos(*it);
	}
	res /= _particles[chain][type].size();

	return res;
}

void GrByInsertion::_set_random_orientation(int chain) {
	LR_matrix R = Utils::get_random_rotation_matrix(2 * M_PI);

	typename vector<BaseParticle *>::iterator it;
	for(it = _particles[chain][0].begin(); it != _particles[chain][0].end(); it++) {
		(*it)->pos = R * _config_info->box->get_abs_pos(*it);
	}
	// we have to move the whole chain
	for(it = _particles[chain][1].begin(); it != _particles[chain][1].end(); it++) {
		(*it)->pos = R * _config_info->box->get_abs_pos(*it);
	}
}

void GrByInsertion::_put_randomly_at_r(int chain, int type, LR_vector &r0, number distance) {
	_set_random_orientation(chain);
	LR_vector com = _get_com(chain, type);
	LR_vector new_com = r0 + Utils::get_random_vector() * distance;
	LR_vector diff = new_com - com;

	typename vector<BaseParticle *>::iterator it;
	for(it = _particles[chain][type].begin(); it != _particles[chain][type].end(); it++) {
		(*it)->pos += diff;
	}
	// we have to move the whole chain
	for(it = _particles[chain][!type].begin(); it != _particles[chain][!type].end(); it++) {
		(*it)->pos += diff;
	}
}

std::string GrByInsertion::get_output_string(llint curr_step) {
	_n_conf++;

	LR_vector coms[2][2];
	for(int c = 0; c < 2; c++)
		for(int t = 0; t < 2; t++)
			coms[c][t] = _get_com(c, t);

	number ref_energy = _config_info->interaction->get_system_energy(_config_info->particles(), _config_info->lists);
	TSPInteraction *_TSP_inter = (TSPInteraction *) _config_info->interaction;
	_TSP_inter->set_only_intra(false);
	for(int i = 0; i < _n_bins; i++) {
		number x0 = i * _bin + _min;
		for(int j = 0; j < _insertions; j++) {
			// random distance between x0 and x0+bin
			number distance = pow(x0 * x0 * x0 + drand48() * (_bin * _bin * _bin + 3. * SQR(_bin) * x0 + 3. * _bin * SQR(x0)), 1. / 3.);
			_put_randomly_at_r(1, _movable, coms[0][_fixed], distance);
			number energy = _TSP_inter->get_system_energy(_config_info->particles(), _config_info->lists);
			// arbitrary threshold
			if(energy < 10000.) {
				number delta_E = energy - ref_energy;
				_inter_hist[_type][i] += exp(-delta_E / _T);
			}
		}
	}
	_TSP_inter->set_only_intra(true);

	stringstream ret;
	number norm = _n_conf * _insertions;
	for(int i = 0; i < _n_bins; i++) {
		number x = (i + 0.5) * _bin + _min;
		number tot_norm = norm;
		ret << x << " ";
		ret << _inter_hist[_type][i] / tot_norm << endl;
	}

	return ret.str();
}
