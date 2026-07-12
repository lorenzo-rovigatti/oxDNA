/*
 * ParticleStress.cpp
 *
 *  Created on: 25/jun/2026
 *      Author: lorenzo
 */
#include "ParticleStress.h"

#include "../Utilities/Utils.h"
#include "../Utilities/oxDNAException.h"

ParticleStress::ParticleStress() :
				BaseObservable() {

}

ParticleStress::~ParticleStress() {

}

void ParticleStress::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputBool(&my_inp, "print_coordinates", &_also_coordinates, 0);
}

std::vector<StressTensor> ParticleStress::stress_tensors() {
	if(!_config_info->interaction->has_custom_stress_tensor()) {
		_config_info->interaction->compute_standard_stress_tensor();
	}

	if(!_config_info->interaction->has_particle_stress_tensor()) {
		throw oxDNAException("ParticleStress requires per-particle stress tensors, "
			"but the current interaction only provides a reduced stress tensor. "
			"If you are running on CUDA, set CUDA_update_particle_stress_tensor = true and ensure that "
			"CUDA_update_stress_tensor_every is set to a value that divides the observable's output interval."
		);
	}

	return _config_info->interaction->particle_stress_tensors();
}

std::string ParticleStress::get_output_string(llint curr_step) {
	auto particle_stress_tensors = stress_tensors();

	std::string to_ret;
	for(int i = 0; i < _config_info->N(); i++) {
		BaseParticle *p = _config_info->particles()[i];
		StressTensor &stress_tensor = particle_stress_tensors[p->index];

		to_ret += Utils::sformat("% .8e % .8e % .8e % .8e % .8e % .8e", stress_tensor[0], stress_tensor[1], stress_tensor[2], stress_tensor[3], stress_tensor[4], stress_tensor[5]);
		if(_also_coordinates) {
			to_ret += Utils::sformat(" % .8e % .8e % .8e", p->pos.x, p->pos.y, p->pos.z);
		}
		to_ret += "\n";
	}

	return to_ret;
}
