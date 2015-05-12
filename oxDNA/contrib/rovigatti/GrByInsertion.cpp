/*
 * GrByInsertion.cpp
 *
 *  Created on: Oct 30, 2013
 *      Author: Lorenzo
 */

#include <sstream>

#include "GrByInsertion.h"
#include "../Utilities/Utils.h"
#include "TSPInteraction.h"

using namespace std;

template<typename number>
GrByInsertion<number>::GrByInsertion() {
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

template<typename number>
GrByInsertion<number>::~GrByInsertion() {
	if(_inter_hist[AA] != NULL) {
		delete[] _inter_hist[AA];
		delete[] _inter_hist[AB];
		delete[] _inter_hist[BB];
	}
}

template<typename number>
void GrByInsertion<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
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

}

template<typename number>
void GrByInsertion<number>::init(ConfigInfo<number> &config_info) {
	BaseObservable<number>::init(config_info);

	if(_max == 0.) _max = *config_info.box_side*0.5;
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

	for(int i = 0; i < *config_info.N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
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

template<typename number>
int GrByInsertion<number>::_get_bin(number sqr_dist) {
	int bin = (sqrt(sqr_dist) - _min) / _bin;
	if(bin >= _n_bins) bin = -1;

	return bin;
}

template<typename number>
LR_vector<number> GrByInsertion<number>::_get_com(int chain, int type) {
	LR_vector<number> res(0., 0., 0.);

	typename vector<BaseParticle<number> *>::iterator it;
	for(it = _particles[chain][type].begin(); it != _particles[chain][type].end(); it++) {
		res += (*it)->get_abs_pos(*this->_config_info.box_side);
	}
	res /= _particles[chain][type].size();

	return res;
}

template<typename number>
void GrByInsertion<number>::_set_random_orientation(int chain) {
	LR_matrix<number> R = Utils::get_random_rotation_matrix<number>(2*M_PI);

	typename vector<BaseParticle<number> *>::iterator it;
	for(it = _particles[chain][0].begin(); it != _particles[chain][0].end(); it++) {
		(*it)->pos = R*(*it)->get_abs_pos(*this->_config_info.box_side);
	}
	// we have to move the whole chain
	for(it = _particles[chain][1].begin(); it != _particles[chain][1].end(); it++) {
		(*it)->pos = R*(*it)->get_abs_pos(*this->_config_info.box_side);
	}
}

template<typename number>
void GrByInsertion<number>::_put_randomly_at_r(int chain, int type, LR_vector<number> &r0, number distance) {
	_set_random_orientation(chain);
	LR_vector<number> com = _get_com(chain, type);
	LR_vector<number> new_com = r0 + Utils::get_random_vector<number>()*distance;
	LR_vector<number> diff = new_com - com;

	typename vector<BaseParticle<number> *>::iterator it;
	for(it = _particles[chain][type].begin(); it != _particles[chain][type].end(); it++) {
		(*it)->pos += diff;
	}
	// we have to move the whole chain
	for(it = _particles[chain][!type].begin(); it != _particles[chain][!type].end(); it++) {
		(*it)->pos += diff;
	}
}

template<typename number>
std::string GrByInsertion<number>::get_output_string(llint curr_step) {
	int N = *this->_config_info.N;
	//number L = *this->_config_info.box_side;
	_n_conf++;

	LR_vector<number> coms[2][2];
	for(int c = 0; c < 2; c++) for(int t = 0; t < 2; t++) coms[c][t] = _get_com(c, t);

	number ref_energy = this->_config_info.interaction->get_system_energy(this->_config_info.particles, N, this->_config_info.lists);
	TSPInteraction<number> *_TSP_inter = (TSPInteraction<number> *) this->_config_info.interaction;
	_TSP_inter->set_only_intra(false);
	for(int i = 0; i < _n_bins; i++) {
		number x0 = i*_bin + _min;
		for(int j = 0; j < _insertions; j++) {
			// random distance between x0 and x0+bin
			number distance = pow(x0*x0*x0 + drand48()*(_bin*_bin*_bin + 3.*SQR(_bin)*x0 + 3.*_bin*SQR(x0)), 1. / 3.);
			_put_randomly_at_r(1, _movable, coms[0][_fixed], distance);
			number energy = _TSP_inter->get_system_energy(this->_config_info.particles, N, this->_config_info.lists);
			// arbitrary threshold
			if(energy < 10000.) {
				number delta_E = energy - ref_energy;
				_inter_hist[_type][i] += exp(-delta_E / _T);
			}
		}
	}
	_TSP_inter->set_only_intra(true);

/*	number factor = 1.;

	// inter
	int bin = _get_bin(coms[0][P_A].sqr_min_image_distance(coms[1][P_A], L));
	if(bin != -1) {
		_inter_hist[AA][bin] += factor;
		_inter_norm[AA] += factor;
	}
	bin = _get_bin(coms[0][P_A].sqr_min_image_distance(coms[1][P_B], L));
	if(bin != -1) {
		_inter_hist[AB][bin] += factor*0.5;
		_inter_norm[AB] += factor*0.5;
	}
	bin = _get_bin(coms[0][P_B].sqr_min_image_distance(coms[1][P_A], L));
	if(bin != -1) {
		_inter_hist[AB][bin] += factor*0.5;
		_inter_norm[AB] += factor*0.5;
	}
	bin = _get_bin(coms[0][P_B].sqr_min_image_distance(coms[1][P_B], L));
	if(bin != -1) {
		_inter_hist[BB][bin] += factor;
		_inter_norm[BB] += factor;
	}*/

	stringstream ret;
	ret << "# r g_AA g_AB g_BB" << endl;
	number norm = _n_conf*_insertions;
	for(int i = 0; i < _n_bins; i++) {
		number x = (i+0.5)*_bin + _min;
		number tot_norm = norm;
		ret << x << " ";
		ret << _inter_hist[_type][i]/tot_norm << endl;
	}

	return ret.str();
}

template class GrByInsertion<float>;
template class GrByInsertion<double>;
