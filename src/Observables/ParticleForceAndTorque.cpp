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

	getInputString(&my_inp, "particle", _particle_string, 1);
	getInputBool(&my_inp, "lab_frame", &_lab_frame, 0);
}

void ParticleForceAndTorque::init() {
	BaseObservable::init();

	auto indexes = Utils::get_particles_from_string(CONFIG_INFO->particles(), _particle_string, "ParticleForceAndTorque");

	for(auto id = indexes.begin(); id != indexes.end(); id++) {
		if(*id < -1 || *id >= _config_info->N()) {
			throw oxDNAException("ParticleForceAndTorque: invalid id %d", *id);
		}

		if(*id == -1) { // all particles
			for(auto p : _config_info->particles()) {
				_particles.insert(p);
			}
		}
		else {
			_particles.insert(CONFIG_INFO->particles()[*id]);
		}
	}
}

std::string ParticleForceAndTorque::get_output_string(llint curr_step) {
	std::string result;
	for(auto p : _particles) {
		LR_vector torque = (_lab_frame) ? p->orientation * p->torque : p->torque;
		result += Utils::sformat("%d %g %g %g %g %g %g\n", p->index, p->force[0], p->force[1], p->force[2], torque[0], torque[1], torque[2]);
	}

	return result;
}
