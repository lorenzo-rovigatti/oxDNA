/*
 * PotentialEnergySelected.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: bonato
 */

#include "PotentialEnergySelected.h"
#include <iostream>

PotentialEnergySelected::PotentialEnergySelected() :
				_split(false) {

}

PotentialEnergySelected::~PotentialEnergySelected() {

}

number PotentialEnergySelected::get_potential_energy() {
	_config_info->interaction->set_is_infinite(false);
	number energy = 0.;
	for(auto n : names) {
		energy += _config_info->interaction->get_system_energy_term(n,_config_info->particles(), _config_info->lists);
	}
	energy /= _config_info->N();

	if(_config_info->interaction->get_is_infinite()) {
		OX_LOG(Logger::LOG_WARNING, "The computation of the total potential energy raised the 'is_infinite' flag");
	}

	return energy;
}

void PotentialEnergySelected::get_settings(input_file &my_inp, input_file &sim_inp) {

	BaseObservable::get_settings(my_inp, sim_inp);	
	getInputBool(&my_inp, "split", &_split, 0);
	bool add = false;
	getInputBool(&my_inp, "back", &add, 0);
	if(add) names.push_back(DNAInteraction::BACKBONE);
	add = false;
	getInputBool(&my_inp, "b_exc", &add, 0);
	if(add) names.push_back(DNAInteraction::BONDED_EXCLUDED_VOLUME);
	add = false;
	getInputBool(&my_inp, "stck", &add, 0);
	if(add) names.push_back(DNAInteraction::STACKING);
	add = false;
	getInputBool(&my_inp, "n_exc", &add, 0);
	if(add) names.push_back(DNAInteraction::NONBONDED_EXCLUDED_VOLUME);
	add = false;
	getInputBool(&my_inp, "hydr", &add, 0);
	if(add) names.push_back(DNAInteraction::HYDROGEN_BONDING);
	add = false;
	getInputBool(&my_inp, "crst", &add, 0);
	if(add) names.push_back(DNAInteraction::CROSS_STACKING);
	add = false;
	getInputBool(&my_inp, "coax", &add, 0);
	if(add) names.push_back(DNAInteraction::COAXIAL_STACKING);
}

std::string PotentialEnergySelected::get_output_string(llint curr_step) {
	if(!_split) {
		number energy = get_potential_energy();

		return Utils::sformat("% 10.6lf", energy);
	}
	else {
		std::string res("");
		
		for(auto n : names) {
			number contrib = _config_info->interaction->get_system_energy_term(n, _config_info->particles(), _config_info->lists);
			contrib = contrib / _config_info->N();
			res = Utils::sformat("%s % 10.6lf", res.c_str(), contrib);
		}

		return res;
	}
}
