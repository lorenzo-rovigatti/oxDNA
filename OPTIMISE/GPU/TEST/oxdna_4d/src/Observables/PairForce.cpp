/*
 * PairForce.cpp
 *
 *  Created on: Jun 27, 2014
 *      Author: majid
 */

#include "PairForce.h"
#include <sstream>
#include <map>

PairForce::PairForce() {

}

PairForce::~PairForce() {

}

void PairForce::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	int tmp = 0;
	_print_all_particles = true;
	if(getInputInt(&my_inp, "particle_id", &tmp, 0) == KEY_FOUND) {
		_print_all_particles = false;
		_particle_id = tmp;
	}
}

std::string PairForce::get_output_string(llint curr_step) {
	BaseParticle *p;
	BaseParticle *q;

	std::stringstream output_str;
	output_str << "#id1 id2 force.x force.y force.z torque.x torque.y torque.z, t = " << curr_step << "\n";

	std::vector<ParticlePair> neighbour_pairs = _config_info->lists->get_potential_interactions();

	for(int i = 0; i < (int) neighbour_pairs.size(); i++) {
		p = neighbour_pairs[i].first;
		q = neighbour_pairs[i].second;

		if(_print_all_particles == false) if(p->index != _particle_id and q->index != _particle_id) continue;

		std::stringstream pair_string;

		// we store and then restore the force since some observables might need the current value	
		LR_vector p_store_f = p->force;
		LR_vector q_store_f = q->force;
		LR_vector p_store_t = p->torque;
		LR_vector q_store_t = q->torque;

		p->force = LR_vector(0.0f, 0.0f, 0.0f);
		q->force = LR_vector(0.0f, 0.0f, 0.0f);
		p->torque = LR_vector(0.0f, 0.0f, 0.0f);
		q->torque = LR_vector(0.0f, 0.0f, 0.0f);

		number pq_interaction = _config_info->interaction->pair_interaction(q, p, true, true);
		if(pq_interaction != (number) 0.f) {
			pair_string << p->index << " " << q->index;
			pair_string << " " << p->force.x << " " << p->force.y << " " << p->force.z << " " << p->torque.x << " " << p->torque.y << " " << p->torque.z << "\n";
			pair_string << q->index << " " << p->index;
			pair_string << " " << q->force.x << " " << q->force.y << " " << q->force.z << " " << q->torque.x << " " << q->torque.y << " " << q->torque.z << "\n";
			output_str << pair_string.str();
		}

		p->force = p_store_f;
		p->torque = p_store_t;
		q->force = q_store_f;
		q->torque = q_store_t;
	}

	return output_str.str();
}
