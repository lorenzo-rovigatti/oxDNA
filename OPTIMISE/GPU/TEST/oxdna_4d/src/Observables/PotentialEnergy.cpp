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
	_config_info->interaction->set_is_infinite(false);
	number energy = _config_info->interaction->get_system_energy(_config_info->particles(), _config_info->lists);
	energy /= _config_info->N();

	if(_config_info->interaction->get_is_infinite()) {
		OX_LOG(Logger::LOG_WARNING, "The computation of the total potential energy raised the 'is_infinite' flag");
	}

	return energy;
}

void PotentialEnergy::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputBool(&my_inp, "split", &_split, 0);
}

std::string PotentialEnergy::get_output_string(llint curr_step) {
	if(!_split) {
		number energy = get_potential_energy();

		return Utils::sformat("% 10.6lf", energy);
	}
	else {
		std::string res("");
		auto energies = _config_info->interaction->get_system_energy_split(_config_info->particles(), _config_info->lists);
		for(auto energy_item : energies) {
			number contrib = energy_item.second / _config_info->N();
			res = Utils::sformat("%s % 10.6lf", res.c_str(), contrib);
		}

		return res;
	}
}
