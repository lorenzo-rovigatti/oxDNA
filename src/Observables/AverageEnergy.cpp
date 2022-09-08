/*
 * AverageEnergy.cpp
 *
 *  Created on: Aug 7, 2018
 *      Author: Erik 
 */

#include "AverageEnergy.h"

#include <sstream>
#include <fstream>

AverageEnergy::AverageEnergy() {

}

AverageEnergy::~AverageEnergy() {

}

void AverageEnergy::init() {
	BaseObservable::init();
	ifstream list;
	list.open(_list_file);
	int n;
	if(!list) {
		OX_LOG(Logger::LOG_ERROR, "average_energy - Can't read particle list file");
		exit(1);
	}
	while(list >> n) {
		_particle_ids.insert(n);
	}
	list.close();
}

void AverageEnergy::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	//Try to load parameters from specific op_file key
	if(getInputString(&my_inp, "nucleotide_list", _list_file, 0) == KEY_FOUND) {
		OX_LOG(Logger::LOG_INFO, "average_energy - loading from particle list file");
	}
	else {
		OX_LOG(Logger::LOG_ERROR, "ERROR: average_energy - No particle list found, aborting");
		exit(1);
	}
}

//sum all energies between particles in the input list and output an average

std::string AverageEnergy::get_output_string(llint curr_step) {

	std::stringstream outstr;

	std::map<int, number> split_energies = _config_info->interaction->get_system_energy_split(_config_info->particles(), _config_info->lists);
	std::vector<ParticlePair> neighbour_pairs = _config_info->lists->get_potential_interactions();

	BaseParticle *p;
	BaseParticle *q;

	number total_energy = 0;
	for(int i = 0; i < (int) neighbour_pairs.size(); i++) {
		p = neighbour_pairs[i].first;
		q = neighbour_pairs[i].second;
		if(std::find(_particle_ids.begin(), _particle_ids.end(), p->get_index()) != _particle_ids.end() and std::find(_particle_ids.begin(), _particle_ids.end(), q->get_index()) != _particle_ids.end()) {
			typename std::map<int, number>::iterator it = split_energies.begin();
			for(; it != split_energies.end(); it++) {
				int k = it->first;
				number k_th_interaction = _config_info->interaction->pair_interaction_term(k, q, p);
				total_energy += k_th_interaction;
			}
		}
	}
	total_energy /= _particle_ids.size();
	outstr << curr_step << ' ' << total_energy << "\n";
	return outstr.str();
}
