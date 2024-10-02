/*
 * BaseList.cpp
 *
 *  Created on: 14 ott 2019
 *      Author: lorenzo
 */

#include "BaseList.h"

void BaseList::get_settings(input_file &inp) {
	char sim_type[512] = "MD";
	getInputString(&inp, "sim_type", sim_type, 0);

	// this chain of ifs is here so that if new sim_types are implemented
	// lists have to know how to manage the associated neighbouring lists
	if(strncmp("MD", sim_type, 512) == 0) _is_MC = false;
	else if(strncmp("MC", sim_type, 512) == 0) _is_MC = true;
	else if(strncmp("MC2", sim_type, 512) == 0) _is_MC = true;
	else if(strncmp("VMMC", sim_type, 512) == 0) _is_MC = true;
	else if(strncmp("PT_VMMC", sim_type, 512) == 0) _is_MC = true;
	else if(strncmp("FFS_MD", sim_type, 512) == 0) _is_MC = false;
	else if(strncmp("min", sim_type, 512) == 0) _is_MC = false;
	else if(strncmp("FIRE", sim_type, 512) == 0) _is_MC = false;
	else throw oxDNAException("BaseList does not know how to handle a '%s' sim_type\n", sim_type);
}

std::vector<BaseParticle *> BaseList::get_all_neighbours(BaseParticle *p) {
	std::vector<BaseParticle *> neighs = get_complete_neigh_list(p);

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
	return neighs;
}

std::vector<ParticlePair > BaseList::get_potential_interactions() {
	// the static here makes sure that the memory allocated every time the function is called is not free'd
	// so that only push_back's that make the vector go beyond their last size trigger re-allocation
	static std::vector<ParticlePair > list;
	list.clear();

	for(auto p : _particles) {
		for(auto q : get_all_neighbours(p)) {
			if(p->index > q->index) {
				list.push_back(ParticlePair(p, q));
			}
		}
	}

	return list;
}
