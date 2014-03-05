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
    _i = -1;
    _j = -1;
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
    if (_i >= *this->_config_info.N || _i < 0) throw oxDNAException ("%s: invalid particle_1 index %d. Aborting", __FILE__, _i);
    if (_j >= *this->_config_info.N || _j < 0) throw oxDNAException ("%s: invalid particle_2 index %d. Aborting", __FILE__, _j);
}

template<typename number>
std::string Distance<number>::get_output_string(llint curr_step) {
    LR_vector<number> dist;
    
    BaseParticle<number> *p = this->_config_info.particles[_i];
    BaseParticle<number> *q = this->_config_info.particles[_j];
    number box = *this->_config_info.box_side;

    if (_PBC) dist = q->pos.minimum_image(p->pos, box);
    else dist = q->get_abs_pos(box) - p->get_abs_pos(box);

    if (_have_dir) dist = _dir * (dist * _dir);

    return Utils::sformat("%14.4lf", sqrt (dist * dist));
}

template<typename number>
void Distance<number>::get_settings (input_file &my_inp, input_file &sim_inp) {
    char tmpstr[256];
    int tmpi;
    
    // particle 1
    getInputInt (&my_inp, "particle_1", &_i, 1); 
    // particle 2
    getInputInt (&my_inp, "particle_2", &_j, 1); 
    
    // optional direction argument
    if (getInputString(&my_inp, "dir", tmpstr, 0) == KEY_FOUND) {
	double x, y, z;
        sscanf(tmpstr, "%lf,%lf,%lf", &x, &y, &z);
        _dir = LR_vector<number> ((number)x, (number)y, (number)z);
        _have_dir = true;
    }
    _dir.normalize();
    
    if (getInputBoolAsInt (&my_inp, "PBC", &tmpi, 0) == KEY_FOUND) _PBC = (tmpi != 0);
    
}

template class Distance<float>;
template class Distance<double>;
