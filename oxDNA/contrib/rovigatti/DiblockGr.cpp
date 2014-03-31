/*
 * DiblockGr.cpp
 *
 *  Created on: Oct 30, 2013
 *      Author: Lorenzo
 */

#include <sstream>

#include "DiblockGr.h"
#include "../Utilities/Utils.h"

template<typename number>
DiblockGr<number>::DiblockGr() {
	_max_dist = 0.;
	_n_bins = 0;
	_n_conf = 0;
	_bin = 0.1;
	_inter_hist[AA] = _inter_hist[AB] = _inter_hist[BB] = _intra_hist = NULL;
}

template<typename number>
DiblockGr<number>::~DiblockGr() {
	if(_intra_hist != NULL) {
		delete[] _inter_hist[AA];
		delete[] _inter_hist[AB];
		delete[] _inter_hist[BB];
		delete[] _intra_hist;
	}
}

template<typename number>
void DiblockGr<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	float tmp;
	if(getInputFloat(&my_inp, "bin", &tmp, 0) == KEY_FOUND) _bin = tmp;
}

template<typename number>
void DiblockGr<number>::init(ConfigInfo<number> &config_info) {
	BaseObservable<number>::init(config_info);

	_max_dist = *config_info.box_side*0.5;
	_n_bins = _max_dist / _bin;
	// rounding
	_bin = _max_dist / _n_bins;

	_inter_hist[AA] = new llint[_n_bins];
	_inter_hist[AB] = new llint[_n_bins];
	_inter_hist[BB] = new llint[_n_bins];
	_intra_hist = new llint[_n_bins];

	for(int i = 0; i < _n_bins; i++) {
		_inter_hist[AA][i] = 0;
		_inter_hist[AB][i] = 0;
		_inter_hist[BB][i] = 0;
		_intra_hist[i] = 0;
	}
}

template<typename number>
int DiblockGr<number>::_get_bin(number sqr_dist) {
	int bin = sqrt(sqr_dist) / _bin;
	if(bin >= _n_bins) return -1;

	return bin;
}

template<typename number>
std::string DiblockGr<number>::get_output_string(llint curr_step) {
	int N = *this->_config_info.N;
	number L = *this->_config_info.box_side;

	LR_vector<number> coms[2][2];
	int counters[2][2] = {0, 0};

	for(int i = 0; i < N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		if(p->strand_id > 1) throw oxDNAException("The system should contain just two chains");
		if(p->type != P_A && p->type != P_B) throw oxDNAException("Only A and B particle types allowed");
		coms[p->strand_id][p->type] += p->get_abs_pos(L);
		counters[p->strand_id][p->type]++;
	}
	for(int i = 0; i < 2; i++) for(int j = 0; j < 2; j++) coms[i][j] /= counters[i][j];

	// inter
	int bin = _get_bin(coms[0][P_A].sqr_min_image_distance(coms[1][P_A], L));
	if(bin != -1) _inter_hist[AA][bin]++;
	bin = _get_bin(coms[0][P_A].sqr_min_image_distance(coms[1][P_B], L));
	if(bin != -1) _inter_hist[AB][bin]++;
	bin = _get_bin(coms[0][P_B].sqr_min_image_distance(coms[1][P_A], L));
	if(bin != -1) _inter_hist[AB][bin]++;
	bin = _get_bin(coms[0][P_B].sqr_min_image_distance(coms[1][P_B], L));
	if(bin != -1) _inter_hist[BB][bin]++;

	// intra
	bin = _get_bin(coms[0][P_A].sqr_min_image_distance(coms[0][P_B], L));
	if(bin != -1) _intra_hist[bin]++;
	bin = _get_bin(coms[1][P_A].sqr_min_image_distance(coms[1][P_B], L));
	if(bin != -1) _intra_hist[bin]++;

	stringstream ret;
	ret << "# r g_AA g_AB g_BB g_intra" << endl;
	_n_conf++;
	number norm = 2.*L*L*L/(2.*4.*M_PI*_bin*_n_conf);
	for(int i = 0; i < _n_bins; i++) {
		number x = (i+0.5)*_bin;
		number tot_norm = norm / SQR(x);
		ret << x << " ";
		ret << _inter_hist[AA][i]*tot_norm << " ";
		ret << _inter_hist[AB][i]*tot_norm/2. << " ";
		ret << _inter_hist[BB][i]*tot_norm << " ";
		ret << _intra_hist[i]*tot_norm/2. << endl;
	}

	return ret.str();
}

template class DiblockGr<float>;
template class DiblockGr<double>;
