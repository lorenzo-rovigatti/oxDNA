/*
 * PairForce.cpp
 *
 *  Created on: Jun 27, 2014
 *      Author: majid
 */

#include "PairForce.h"
#include <sstream>
#include <map>

template<typename number>
PairForce<number>::PairForce() {

}

template<typename number>
PairForce<number>::~PairForce() {

}

template<typename number>
void PairForce<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	int tmp = 0;
	_print_all_particles = true;
	if (getInputInt(&my_inp,"particle_id", &tmp, 0) == KEY_FOUND){
		_print_all_particles = false;
		_particle_id = tmp;
	}
}

template<typename number>
std::string PairForce<number>::get_output_string(llint curr_step) {
	BaseParticle<number> *p;
	BaseParticle<number> *q;

	std::stringstream output_str;
	output_str << "#id1 id2 force.x force.y force.z torque.x torque.y torque.z, t = " << curr_step << "\n";

	std::vector<ParticlePair<number> > neighbour_pairs = this->_config_info.lists->get_potential_interactions();

	for (int i = 0; i < (int)neighbour_pairs.size(); i++) {
		p = neighbour_pairs[i].first;
		q = neighbour_pairs[i].second;
	
		if (_print_all_particles == false)
			if (p->index != _particle_id and q->index != _particle_id)
				continue;

		std::stringstream pair_string;
		
		// we store and then restore the force since some observables might need the current value	
		LR_vector<number> p_store_f = p->force;
		LR_vector<number> q_store_f = q->force;
		LR_vector<number> p_store_t = p->torque;
		LR_vector<number> q_store_t = q->torque;

		p->force = LR_vector<number> (0.0f, 0.0f, 0.0f);
		q->force = LR_vector<number> (0.0f, 0.0f, 0.0f);
		p->torque = LR_vector<number> (0.0f, 0.0f, 0.0f);
		q->torque = LR_vector<number> (0.0f, 0.0f, 0.0f);

		number pq_interaction = this->_config_info.interaction->pair_interaction (q, p, NULL, true);
		if(pq_interaction != (number) 0.f) {
			pair_string << p->index << " " << q->index;
			pair_string << " " <<  p->force.x << " " <<  p->force.y << " " <<  p->force.z << " "  <<  p->torque.x << " " <<  p->torque.y  << " " <<  p->torque.z << "\n";
			pair_string << q->index << " " << p->index;
			pair_string << " " <<  q->force.x << " " <<  q->force.y << " " <<  q->force.z << " "  <<  q->torque.x << " " <<  q->torque.y  << " " <<  q->torque.z << "\n";
			output_str << pair_string.str();
		}
		
		p->force = p_store_f;
		p->torque = p_store_t;
		q->force = q_store_f;
		q->torque = q_store_t;
	}

    return output_str.str();
}

template class PairForce<float>;
template class PairForce<double>;

