/*
 * GenericGrByInsertion.cpp
 *
 *  Created on: 12/jan/2017
 *      Author: lorenzo
 */
#include <sstream>

#include "GenericGrByInsertion.h"
#include "Utilities/Utils.h"

using namespace std;

template<typename number>
GenericGrByInsertion<number>::GenericGrByInsertion() {
	_n_bins = 0;
	_n_conf = 0;
	_bin = 0.1;
	_T = 0.;
	_gr = NULL;
	_insertions = 100;
	_min = 0.0;
	_max = 15;
	_second_starts_from = -1;
}

template<typename number>
GenericGrByInsertion<number>::~GenericGrByInsertion() {
	if(_gr != NULL) delete[] _gr;
}

template<typename number>
void GenericGrByInsertion<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	float tmp;
	if(getInputFloat(&my_inp, "bin", &tmp, 0) == KEY_FOUND) _bin = tmp;
	
	getInputInt(&my_inp, "insertions", &_insertions, 0);

	char tmp_T[256];
	getInputString(&sim_inp, "T", tmp_T, 1);
	_T = Utils::get_temperature<number>(tmp_T);
	
	if(getInputFloat(&my_inp, "gr_min", &tmp, 0) == KEY_FOUND) _min = tmp;
	if(getInputFloat(&my_inp, "gr_max", &tmp, 0) == KEY_FOUND) _max = tmp;
	getInputInt(&my_inp, "second_starts_from", &_second_starts_from, 0);
	
	CHECK_BOX("GrByInsertion", my_inp);
}

template<typename number>
int GenericGrByInsertion<number>::_get_bin(number sqr_dist) {
	int bin = (sqrt(sqr_dist) - _min) / _bin;
	if(bin >= _n_bins) bin = -1;
	
	return bin;
}

template<typename number>
void GenericGrByInsertion<number>::init(ConfigInfo<number> &config_info) {
	BaseObservable<number>::init(config_info);
	int N = *config_info.N;
	
	if(_max == 0.) _max = config_info.box->box_sides().x*0.5;
	_n_bins = (_max - _min) / _bin;
	
	_gr = new number[_n_bins];
	
	for(int i = 0; i < _n_bins; i++) _gr[i] = 0.;
	
	for(int i = 0; i < N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		if(_second_starts_from == -1) _particles[(int)(i/(N/2))].push_back(p);
		else {
			int idx = 1;
			if(i < _second_starts_from) idx = 0;
			_particles[idx].push_back(p);
		}
	}
}

template<typename number>
LR_vector<number> GenericGrByInsertion<number>::_get_com(int obj) {
	LR_vector<number> res(0., 0., 0.);
	
	typename vector<BaseParticle<number> *>::iterator it;
	for(it = _particles[obj].begin(); it != _particles[obj].end(); it++) {
		res += this->_config_info.box->get_abs_pos(*it);
	}
	res /= _particles[obj].size();

	return res;
}

template<typename number>
void GenericGrByInsertion<number>::_set_random_orientation(int obj) {
	LR_matrix<number> R = Utils::get_random_rotation_matrix<number>(2*M_PI);
	
	typename vector<BaseParticle<number> *>::iterator it;
	
	for(it = _particles[obj].begin(); it != _particles[obj].end(); it++) {
		(*it)->pos = R*this->_config_info.box->get_abs_pos(*it);
		(*it)->orientation = R*(*it)->orientation;
		(*it)->set_positions();
		(*it)->orientationT = (*it)->orientation.get_transpose();
	}
}

template<typename number>
void GenericGrByInsertion<number>::_put_randomly_at_r(int obj, LR_vector<number> &r0, number distance) {
	_set_random_orientation(obj);
	LR_vector<number> com = _get_com(obj);
	LR_vector<number> new_com = r0 + Utils::get_random_vector<number>()*distance;
	LR_vector<number> diff = new_com - com;
	
	typename vector<BaseParticle<number> *>::iterator it;
	
	for(it = _particles[obj].begin(); it != _particles[obj].end(); it++) (*it)->pos += diff;
}

template<typename number>
std::string GenericGrByInsertion<number>::get_output_string(llint step) {
	int N = *this->_config_info.N;

	_n_conf++;
	
	LR_vector<number> coms[2];
	for(int t = 0; t < 2; t++) coms[t] = _get_com(t);
		
	number ref_energy = this->_config_info.interaction->get_system_energy(this->_config_info.particles, N, this->_config_info.lists);
	for(int i = 0; i < _n_bins; i++) {
		number x0 = i*_bin + _min;
		for(int j = 0; j < _insertions; j++) {
			
			_put_randomly_at_r(0, coms[1], x0);
			this->_config_info.lists->global_update(true);
			number energy = this->_config_info.interaction->get_system_energy(this->_config_info.particles, N, this->_config_info.lists);
			// arbitrary threshold
			if(energy / N < 1000.) {
				number delta_E = energy - ref_energy;
				number boltzmann_factor = exp(-delta_E / _T);
				
				// center center g(r)
				_gr[i] += boltzmann_factor;				
			}
		}
	}
	// now we move the objects very far apart
	_put_randomly_at_r (0, coms[1], 100);
	this->_config_info.lists->global_update(true);
	
	stringstream myret;

	number norm = _n_conf * _insertions;
	for(int i = 0; i < _n_bins; i++) {
		number x = i*_bin + _min;
		myret << x << " " << _gr[i] / norm << endl;
	}

	return myret.str();
}

template class GenericGrByInsertion<float>;
template class GenericGrByInsertion<double>;
