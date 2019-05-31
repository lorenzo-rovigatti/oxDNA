/*
 * PotentialEnergy.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#include "PotentialEnergy.h"

template<typename number>
PotentialEnergy<number>::PotentialEnergy(): _split(false) {

}

template<typename number>
PotentialEnergy<number>::~PotentialEnergy() {

}

template<typename number>
number PotentialEnergy<number>::get_potential_energy() {
	number energy = this->_config_info.interaction->get_system_energy(this->_config_info.particles, *this->_config_info.N, this->_config_info.lists);
	energy /= *this->_config_info.N;

	return energy;
}

template<typename number>
void PotentialEnergy<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	getInputBool(&my_inp, "split", &_split, 0);
}

template<typename number>
std::string PotentialEnergy<number>::get_output_string(llint curr_step) {
	if(!_split) {
		number energy = get_potential_energy();

		return Utils::sformat("% 10.6lf", energy);
	}
	else {
		string res("");
		map<int, number> energies = this->_config_info.interaction->get_system_energy_split(this->_config_info.particles, *this->_config_info.N, this->_config_info.lists);
		for(typename map<int, number>::iterator it = energies.begin(); it != energies.end(); it++) {
			number contrib = it->second / *this->_config_info.N;
			res = Utils::sformat("%s % 10.6lf", res.c_str(), contrib);
		}

		return res;
	}
}

template class PotentialEnergy<float>;
template class PotentialEnergy<double>;
