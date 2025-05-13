/*
 * ParticleForceAndTorque.cpp
 *
 *  Created on: May 13, 2025
 *      Author: lorenzo
 */

#include "ParticleForceAndTorque.h"

ParticleForceAndTorque::ParticleForceAndTorque() {

}

ParticleForceAndTorque::~ParticleForceAndTorque() {

}

void ParticleForceAndTorque::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	std::string p_list;
	getInputString(&my_inp, "particle", p_list, 1);

	_indexes = Utils::get_particles_from_string(CONFIG_INFO->particles(), p_list, "ParticleForceAndTorque");
}

void ParticleForceAndTorque::init() {
	BaseObservable::init();

	for(auto id = _indexes.begin(); id != _indexes.end(); id++) {
		if(*id < 0 || *id >= _config_info->N()) {
			throw oxDNAException("ParticleForceAndTorque: invalid id %d", *id);
		}
		_particles.push_back(CONFIG_INFO->particles()[*id]);
	}
}

std::string ParticleForceAndTorque::get_output_string(llint curr_step) {
	std::string result;
	for(auto p : _particles) {
		LR_vector torque = p->orientation * p->torque;
		result += Utils::sformat("%d %g %g %g %g %g %g\n", p->index, p->force[0], p->force[1], p->force[2], torque[0], torque[1], torque[2]);
	}

	return result;
}
