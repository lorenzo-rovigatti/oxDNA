/*
 * meta_utils.cpp
 *
 *  Created on: Apr 3, 2022
 *      Author: lorenzo
 */

#include "meta_utils.h"

namespace meta {

LR_vector particle_list_com(const std::vector<BaseParticle *> &list, BaseBox *box_ptr) {
	auto sum_abs_pos = [&](LR_vector sum, BaseParticle *p) { return sum + box_ptr->get_abs_pos(p); };
	return std::accumulate(list.begin(), list.end(), LR_vector(), sum_abs_pos) / (number) list.size();
}

std::tuple<std::vector<int>, std::vector<BaseParticle *>> get_particle_lists(input_file &inp, std::string key, std::vector<BaseParticle *> &particles, std::string description) {
	std::string p_string;
	getInputString(&inp, key.c_str(), p_string, 1);

	std::vector<int> idx_vector;
	std::vector<BaseParticle *> p_vector;
	idx_vector = Utils::get_particles_from_string(particles, p_string, description);
	for(auto p_idx : idx_vector) {
		p_vector.push_back(particles[p_idx]);
	}

	return std::make_tuple(idx_vector, p_vector);
}

}
