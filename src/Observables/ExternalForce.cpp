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

	getInputString(&my_inp, "particles", _particle_ids, 0);
}

void ExternalForce::init() {
	std::vector<int> ids = Utils::get_particles_from_string(_config_info->particles(), _particle_ids, "ExternalForce observable");

	if(ids.size() == 0) {
		_particles = _config_info->particles();
	}
	else {
		for(auto p : _config_info->particles()) {
			if(std::find(ids.begin(), ids.end(), p->index) != std::end(ids)) {
				_particles.push_back(p);
			}
		}
	}
	_force_averages.resize(_particles.size());
}

void ExternalForce::update_data(llint curr_step) {
	for(uint i = 0; i < _particles.size(); i++) {
		auto p = _particles[i];
		p->set_initial_forces(curr_step, _config_info->box);
		_force_averages[i] += p->force;
	}
	_times_updated++;
}

std::string ExternalForce::get_output_string(llint curr_step) {
	std::string output = Utils::sformat("# t = %lld\n", curr_step);

	for(auto &force : _force_averages) {
		force /= _times_updated;
		output += Utils::sformat("%lf %lf %lf\n", force.x, force.y, force.z);
		force = LR_vector();
	}
	_times_updated = 0;

	return output;
}
