/*
 * PotentialEnergy.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#include "PotentialEnergy.h"

PotentialEnergy::PotentialEnergy() :
				_split(false) {

}

PotentialEnergy::~PotentialEnergy() {

}

number PotentialEnergy::get_potential_energy() {
	number energy = _config_info->interaction->get_system_energy(_config_info->particles, *_config_info->N, _config_info->lists);
	energy /= *_config_info->N;

	return energy;
}

void PotentialEnergy::get_settings(input_file &my_inp, input_file &sim_inp) {
	getInputBool(&my_inp, "split", &_split, 0);
}

std::string PotentialEnergy::get_output_string(llint curr_step) {
	if(!_split) {
		number energy = get_potential_energy();

		return Utils::sformat("% 10.6lf", energy);
	}
	else {
		std::string res("");
		std::map<int, number> energies = _config_info->interaction->get_system_energy_split(_config_info->particles, *_config_info->N, _config_info->lists);
		for(auto it = energies.begin(); it != energies.end(); it++) {
			number contrib = it->second / *_config_info->N;
			res = Utils::sformat("%s % 10.6lf", res.c_str(), contrib);
		}

		return res;
	}
}
