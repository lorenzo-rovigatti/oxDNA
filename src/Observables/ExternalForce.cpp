/*
 * ExternalForce.cpp
 *
 *  Created on: Feb 2, 2024
 *      Author: rovigatti
 */

#include "ExternalForce.h"

ExternalForce::ExternalForce() {

}

ExternalForce::~ExternalForce() {

}

void ExternalForce::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	std::string ids;
	if(getInputString(&my_inp, "ids", ids, 0) == KEY_FOUND) {
		_ids = Utils::split(ids, ',');
	}
}

void ExternalForce::update_data(llint curr_step) {
	static bool initialised = false;
	if(!initialised) {
		if(_ids.size() == 0) {
			_force_averages.resize(_config_info->forces.size());
			_ext_forces = _config_info->forces;
		}
		else {
			_force_averages.resize(_ids.size());
			for(auto id : _ids) {
				_ext_forces.push_back(_config_info->get_force_by_id(id));
			}
			_ext_forces.resize(_ids.size());
		}
		initialised = true;
	}
}

std::string ExternalForce::get_output_string(llint curr_step) {
	return Utils::sformat("TO BE IMPLEMENTED");
}
