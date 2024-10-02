/*
 * PatchyBonds.cpp
 *
 *  Created on: 16/aug/2020
 *      Author: lorenzo
 */

#include <sstream>
#include <fstream>
#include <iostream>
#include <string>

#include "PatchyBonds.h"
#include "Interactions/PatchyInteraction.h"
#include "../Interactions/FSInteraction.h"

using namespace std;

PatchyBonds::PatchyBonds() :
				Configuration() {

}

PatchyBonds::~PatchyBonds() {

}

void PatchyBonds::get_settings(input_file &my_inp, input_file &sim_inp) {
	Configuration::get_settings(my_inp, sim_inp);
}

std::string PatchyBonds::_headers(llint step) {
	std::stringstream headers;

	headers << "# " << "step " << _config_info->curr_step << " N " << _config_info->N() << endl;

	return headers.str();
}

string PatchyBonds::_configuration(llint step) {
	stringstream conf;
	conf.precision(15);

	_config_info->interaction->get_system_energy(_config_info->particles(), _config_info->lists);

	for(auto p : _config_info->particles()) {
		PatchyParticle *pp = static_cast<PatchyParticle *>(p);

		std::vector<std::vector<int>> bonded_with(pp->N_int_centers());

		for(auto &bond : pp->bonds) {
			bonded_with[bond.p_patch].push_back(bond.other->index);
		}

		for(auto &patch_bonds : bonded_with) {
			conf << patch_bonds.size() << " ";
		}
		conf << std::endl;

		for(auto &patch_bonds : bonded_with) {
			for(auto bonded_id : patch_bonds) {
				conf << (bonded_id + 1) << " ";
			}
		}
		conf << std::endl;
	}

	return conf.str();
}
