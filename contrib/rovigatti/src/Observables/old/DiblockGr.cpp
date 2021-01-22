/*
 * DiblockGr.cpp
 *
 *  Created on: Oct 30, 2013
 *      Author: Lorenzo
 */

#include <sstream>

#include "DiblockGr.h"
#include "Utilities/Utils.h"

DiblockGr::DiblockGr() {
	_max_dist = 0.;
	_n_bins = 0;
	_n_conf = 0;
	_bin = 0.1;
	_biased = false;
	_T = 0.;
	_inter_hist[AA] = _inter_hist[AB] = _inter_hist[BB] = _intra_hist = NULL;
	_inter_norm[AA] = _inter_norm[AB] = _inter_norm[BB] = _intra_norm = 0.;
	_only_intra = 0;
}

DiblockGr::~DiblockGr() {
	if(_intra_hist != NULL) {
		delete[] _inter_hist[AA];
		delete[] _inter_hist[AB];
		delete[] _inter_hist[BB];
		delete[] _intra_hist;
	}
}

void DiblockGr::get_settings(input_file &my_inp, input_file &sim_inp) {
	float tmp;
	if(getInputFloat(&my_inp, "bin", &tmp, 0) == KEY_FOUND) _bin = tmp;
	int biased = 0;
	getInputBoolAsInt(&my_inp, "biased", &biased, 0);
	_biased = (bool) biased;
	int only_intra = 0;
	getInputBoolAsInt(&my_inp, "only_intra", &only_intra, 0);
	_only_intra = (bool) only_intra;

	CHECK_BOX("DiblockGr", my_inp);

	if(_biased) {
		_force_energy.get_settings(my_inp, sim_inp);
		getInputFloat(&sim_inp, "T", &tmp, 1);
		_T = tmp;
		OX_LOG(Logger::LOG_INFO, "Biased simulation, DiblockGr will automatically unbias the g(r)");
	}
}

void DiblockGr::init() {
	BaseObservable::init();

	_max_dist = _config_info->box->box_sides().x * 0.5;
	_n_bins = _max_dist / _bin;
	// rounding
	_bin = _max_dist / _n_bins;

	_inter_hist[AA] = new number[_n_bins];
	_inter_hist[AB] = new number[_n_bins];
	_inter_hist[BB] = new number[_n_bins];
	_intra_hist = new number[_n_bins];

	for(int i = 0; i < _n_bins; i++) {
		_inter_hist[AA][i] = 0.;
		_inter_hist[AB][i] = 0.;
		_inter_hist[BB][i] = 0.;
		_intra_hist[i] = 0.;
	}

	if(_biased) _force_energy.init(config_info);
}

int DiblockGr::_get_bin(number sqr_dist) {
	int bin = sqrt(sqr_dist) / _bin;
	if(bin >= _n_bins) bin = -1;

	return bin;
}

std::string DiblockGr::get_output_string(llint curr_step) {
	int N = _config_info->N();
	_n_conf++;

	LR_vector coms[2][2];
	int counters[2][2] = { { 0, 0 }, { 0, 0 } };

	for(int i = 0; i < N; i++) {
		BaseParticle *p = _config_info->particles()[i];
		if(p->strand_id > 1) throw oxDNAException("The system should contain just two chains");
		if(p->type != P_A && p->type != P_B) throw oxDNAException("Only particles of type A and B are allowed");
		coms[p->strand_id][p->type] += _config_info->box->get_abs_pos(p);
		counters[p->strand_id][p->type]++;
	}
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
			coms[i][j] /= counters[i][j];

	number factor = 2.;
	if(_biased) {
		// we need N because the ForceEnergy observable computes the energy per particle and
		// 0.5 because we have to put two traps to exert a "real", newtonian force between
		// the two particles
		number f_e = atof(_force_energy.get_output_string(curr_step).c_str()) * N * 0.5;
		number val = exp(f_e / _T);
		factor *= val;
	}

	// inter
	int bin;
	if(!_only_intra) {
		bin = _get_bin(_config_info->box->sqr_min_image_distance(coms[1][P_A], coms[0][P_A]));
		if(bin != -1) {
			_inter_hist[AA][bin] += factor;
			_inter_norm[AA] += factor;
		}
		bin = _get_bin(_config_info->box->sqr_min_image_distance(coms[1][P_B], coms[0][P_A]));
		if(bin != -1) {
			_inter_hist[AB][bin] += factor * 0.5;
			_inter_norm[AB] += factor * 0.5;
		}
		bin = _get_bin(_config_info->box->sqr_min_image_distance(coms[1][P_A], coms[0][P_B]));
		if(bin != -1) {
			_inter_hist[AB][bin] += factor * 0.5;
			_inter_norm[AB] += factor * 0.5;
		}
		bin = _get_bin(_config_info->box->sqr_min_image_distance(coms[1][P_B], coms[0][P_B]));
		if(bin != -1) {
			_inter_hist[BB][bin] += factor;
			_inter_norm[BB] += factor;
		}
	}

	// intra
	bin = _get_bin(_config_info->box->sqr_min_image_distance(coms[0][P_B], coms[0][P_A]));
	if(bin != -1) {
		_intra_hist[bin] += factor * 0.5;
		_intra_norm += factor * 0.5;
	}

	if(!_only_intra) {
		bin = _get_bin(_config_info->box->sqr_min_image_distance(coms[1][P_B], coms[1][P_A]));
		if(bin != -1) {
			_intra_hist[bin] += factor * 0.5;
			_intra_norm += factor * 0.5;
		}
	}

	std::stringstream ret;
	ret << "# r g_AA g_AB g_BB g_intra" << std::endl;
	number norm = (4. * M_PI) / (3. * _config_info->box->V()) * 2;
	for(int i = 0; i < _n_bins; i++) {
		number x0 = i * _bin;
		number x1 = (i + 1) * _bin;
		number vb = (x1 * x1 * x1 - x0 * x0 * x0);
		number tot_norm = norm * vb;
		ret << x0 + 0.5 * _bin << " ";
		if(!_only_intra) {
			ret << _inter_hist[AA][i] / (tot_norm * _inter_norm[AA]) << " ";
			ret << _inter_hist[AB][i] / (tot_norm * _inter_norm[AB]) << " ";
			ret << _inter_hist[BB][i] / (tot_norm * _inter_norm[BB]) << " ";
		}
		ret << _intra_hist[i] / _intra_norm << " " << _intra_hist[i] / (_intra_norm * tot_norm) << std::endl;
	}

	return ret.str();
}
