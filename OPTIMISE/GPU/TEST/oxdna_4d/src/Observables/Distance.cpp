/*
 * Distance.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#include "Distance.h"
#include "../Utilities/Utils.h"

Distance::Distance() {
	_PBC = true;
	_have_dir = false;
	_dir = LR_vector(1., 1., 1.);
}

Distance::~Distance() {

}

void Distance::init() {
	BaseObservable::init();

	int N = _config_info->N();
	std::vector<BaseParticle *> &particles = _config_info->particles();

	std::vector<int> p1_indexes = Utils::get_particles_from_string(particles, _p1_string, "Distance observable");
	for(auto idx: p1_indexes) {
		_check_index(idx, N);
		_p1_list.insert(particles[idx]);
	}

	std::vector<int> p2_indexes = Utils::get_particles_from_string(particles, _p2_string, "Distance observable");
	for(auto idx: p2_indexes) {
		_check_index(idx, N);
		_p2_list.insert(particles[idx]);
	}
}

std::string Distance::get_output_string(llint curr_step) {
	LR_vector dist;
	LR_vector p1_com, p2_com;
	for(auto particle : _p1_list) {
		p1_com += _config_info->box->get_abs_pos(particle);
	}
	p1_com /= _p1_list.size();

	for(auto particle : _p2_list) {
		p2_com += _config_info->box->get_abs_pos(particle);
	}
	p2_com /= _p2_list.size();

	if(_PBC) dist = _config_info->box->min_image(p1_com, p2_com);
	else dist = p2_com - p1_com;

	number sign = 1.;
	if(_have_dir) {
		dist = _dir * (dist * _dir);
		sign = copysign(1., dist * _dir);
	}

	return Utils::sformat("%14.4lf", sign * sqrt(dist * dist));
}

void Distance::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	// particle 1
	getInputString(&my_inp, "particle_1", _p1_string, 1);
	// particle 2
	getInputString(&my_inp, "particle_2", _p2_string, 1);

	// optional direction argument
	char tmpstr[256];
	if(getInputString(&my_inp, "dir", tmpstr, 0) == KEY_FOUND) {
		double x, y, z;
		sscanf(tmpstr, "%lf,%lf,%lf", &x, &y, &z);
		_dir = LR_vector((number) x, (number) y, (number) z);
		_have_dir = true;
	}
	_dir.normalize();

	getInputBool(&my_inp, "PBC", &_PBC, 0);

}
