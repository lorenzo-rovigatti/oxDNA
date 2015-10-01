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
    vector<string> spl = Utils::split(_p1_string, ',');
	for(vector<string>::iterator it = spl.begin(); it != spl.end(); it++) {
		int index = atoi(it->c_str());
		_check_index(index, N);
		_p1_list.insert(particles[index]);
	}

	spl = Utils::split(_p2_string, ',');
	for(vector<string>::iterator it = spl.begin(); it != spl.end(); it++) {
		int index = atoi(it->c_str());
		_check_index(index, N);
		_p2_list.insert(particles[index]);
	}
}

template<typename number>
std::string Distance<number>::get_output_string(llint curr_step) {
    LR_vector<number> dist;
    number box = *this->_config_info.box_side;
    
    LR_vector<number> p1_com, p2_com;
	for(typename set<BaseParticle<number> *>::iterator it = _p1_list.begin(); it != _p1_list.end(); it++) {
		p1_com += (*it)->get_abs_pos(box);
	}
	p1_com /= _p1_list.size();

	for(typename set<BaseParticle<number> *>::iterator it = _p2_list.begin(); it != _p2_list.end(); it++) {
		p2_com += (*it)->get_abs_pos(box);
	}
	p2_com /= _p2_list.size();

    if (_PBC) dist = p2_com.minimum_image(p1_com, box);
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
