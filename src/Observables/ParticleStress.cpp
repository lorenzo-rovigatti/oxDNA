/*
 * ParticleStress.cpp
 *
 *  Created on: 25/jun/2026
 *      Author: lorenzo
 */
#include "ParticleStress.h"

#include "../Utilities/Utils.h"

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
	std::vector<StressTensor> stress_tensors(_config_info->N());

	if(_config_info->interaction->has_custom_stress_tensor()) {
		for(int i = 0; i < _config_info->N(); i++) {
			BaseParticle *p = _config_info->particles()[i];
			StressTensor &stress_tensor = stress_tensors[p->index];

			stress_tensor[0] += p->pos.x * p->force.x;
			stress_tensor[1] += p->pos.y * p->force.y;
			stress_tensor[2] += p->pos.z * p->force.z;
			stress_tensor[3] += p->pos.x * p->force.y;
			stress_tensor[4] += p->pos.x * p->force.z;
			stress_tensor[5] += p->pos.y * p->force.z;

			stress_tensor[0] += SQR(p->vel.x);
			stress_tensor[1] += SQR(p->vel.y);
			stress_tensor[2] += SQR(p->vel.z);
			stress_tensor[3] += p->vel.x * p->vel.y;
			stress_tensor[4] += p->vel.x * p->vel.z;
			stress_tensor[5] += p->vel.y * p->vel.z;
		}
	}
	else {
		static std::vector<LR_vector> old_forces, old_torques;
		_config_info->interaction->begin_energy_and_force_computation();

		// pair_interaction will change these vectors, but we still need them in the next
		// first integration step. For this reason we copy and then restore their values
		// after the calculation
		old_forces.resize(_config_info->N());
		old_torques.resize(_config_info->N());
		for(int i = 0; i < _config_info->N(); i++) {
			old_forces[i] = _config_info->particles()[i]->force;
			old_torques[i] = _config_info->particles()[i]->torque;
		}
		
		for(int i = 0; i < _config_info->N(); i++) {
			BaseParticle *p = _config_info->particles()[i];
			StressTensor &stress_tensor = stress_tensors[p->index];
			std::vector<BaseParticle *> neighs = _config_info->lists->get_neigh_list(p);

			std::set<BaseParticle *> bonded_neighs;
			for(auto &pair : p->affected) {
				if(pair.first != p) {
					bonded_neighs.insert(pair.first);
				}
				if(pair.second != p) {
					bonded_neighs.insert(pair.second);
				}
			}

			neighs.insert(neighs.end(), bonded_neighs.begin(), bonded_neighs.end());

			for(auto q : neighs) {
				if(p->index > q->index) {
					LR_vector r = _config_info->box->min_image(p->pos, q->pos);
					_config_info->interaction->set_computed_r(r);

					p->force = LR_vector();
					_config_info->interaction->pair_interaction(p, q, false, true);

					stress_tensor[0] -= r.x * p->force.x;
					stress_tensor[1] -= r.y * p->force.y;
					stress_tensor[2] -= r.z * p->force.z;
					stress_tensor[3] -= r.x * p->force.y;
					stress_tensor[4] -= r.x * p->force.z;
					stress_tensor[5] -= r.y * p->force.z;
				}
			}

			stress_tensor[0] += SQR(p->vel.x);
			stress_tensor[1] += SQR(p->vel.y);
			stress_tensor[2] += SQR(p->vel.z);
			stress_tensor[3] += p->vel.x * p->vel.y;
			stress_tensor[4] += p->vel.x * p->vel.z;
			stress_tensor[5] += p->vel.y * p->vel.z;
		}

		for(int i = 0; i < _config_info->N(); i++) {
			_config_info->particles()[i]->force = old_forces[i];
			_config_info->particles()[i]->torque = old_torques[i];
		}
	}

	return stress_tensors;
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
