/*
 * Distance.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#include "Distance.h"
#include "../Utilities/Utils.h"

template<typename number>
Distance<number>::Distance() {
    _PBC = true;
    _have_dir = false;
    _dir = LR_vector<number> (1., 1., 1.);
}

template<typename number>
Distance<number>::~Distance() {

}

template<typename number>
void Distance<number>::init(ConfigInfo<number> &config_info) {
    BaseObservable<number>::init(config_info);

    int N = *config_info.N;
    BaseParticle<number> **particles = config_info.particles;

    vector<int> p1_indexes = Utils::getParticlesFromString(particles, N, _p1_string, "Distance observable");
    for(vector<int>::iterator it = p1_indexes.begin(); it != p1_indexes.end(); it++) {
    	_check_index(*it, N);
    	_p1_list.insert(particles[*it]);
    }

    vector<int> p2_indexes = Utils::getParticlesFromString(particles, N, _p2_string, "Distance observable");
    for(vector<int>::iterator it = p2_indexes.begin(); it != p2_indexes.end(); it++) {
    	_check_index(*it, N);
    	_p2_list.insert(particles[*it]);
    }
}

template<typename number>
std::string Distance<number>::get_output_string(llint curr_step) {
    LR_vector<number> dist;
    LR_vector<number> p1_com, p2_com;
	for(typename set<BaseParticle<number> *>::iterator it = _p1_list.begin(); it != _p1_list.end(); it++) {
		p1_com += this->_config_info.box->get_abs_pos(*it);
	}
	p1_com /= _p1_list.size();

	for(typename set<BaseParticle<number> *>::iterator it = _p2_list.begin(); it != _p2_list.end(); it++) {
		p2_com += this->_config_info.box->get_abs_pos(*it);
	}
	p2_com /= _p2_list.size();

    if (_PBC) dist = this->_config_info.box->min_image(p1_com,p2_com);
    else dist = p2_com - p1_com;

	number sign  = 1.;
	if (_have_dir) {
		dist = _dir * (dist * _dir);
		sign = copysign (1., dist * _dir);
	}

	return Utils::sformat("%14.4lf", sign*sqrt(dist*dist));
}

template<typename number>
void Distance<number>::get_settings (input_file &my_inp, input_file &sim_inp) {
	// particle 1
	getInputString(&my_inp, "particle_1", _p1_string, 1);
	// particle 2
	getInputString(&my_inp, "particle_2", _p2_string, 1);

    // optional direction argument
    char tmpstr[256];
    if (getInputString(&my_inp, "dir", tmpstr, 0) == KEY_FOUND) {
	double x, y, z;
        sscanf(tmpstr, "%lf,%lf,%lf", &x, &y, &z);
        _dir = LR_vector<number> ((number)x, (number)y, (number)z);
        _have_dir = true;
    }
    _dir.normalize();
    
    getInputBool(&my_inp, "PBC", &_PBC, 0);
    
}

template class Distance<float>;
template class Distance<double>;


