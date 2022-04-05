/*
 * Bonds.cpp
 *
 *  Created on: 19/jul/2021
 *      Author: lorenzo
 */

#include <sstream>
#include <fstream>
#include <iostream>
#include <string>

#include "Bonds.h"
#include "../Interactions/DetailedPolymerSwapInteraction.h"

using namespace std;

Bonds::Bonds() :
				BaseObservable() {

}

Bonds::~Bonds() {

}

void Bonds::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputNumber(&my_inp, "bond_threshold", &_bond_threshold, 0);
	getInputInt(&my_inp, "bond_energy_term_id", &_energy_term_id, 0);
}

void Bonds::init() {
	BaseObservable::init();

	_bonds.resize(_config_info->N());
}

string Bonds::get_output_string(llint step) {
	stringstream outstr;
	outstr << "# step " << step << "\n";

	for(int i = 0; i < _config_info->N(); i++) {
		_bonds[i].clear();
	}

	// compute the bonding pattern
	vector<ParticlePair> inter_pairs = _config_info->lists->get_potential_interactions();
	_config_info->interaction->begin_energy_computation();

	DetailedPolymerSwapInteraction *interaction = dynamic_cast<DetailedPolymerSwapInteraction *>(_config_info->interaction);
	if(interaction != nullptr) {
		interaction->no_three_body = true;
	}

	for(typename vector<ParticlePair>::iterator it = inter_pairs.begin(); it != inter_pairs.end(); it++) {
		number energy;
		if(_energy_term_id == -1) {
			energy = _config_info->interaction->pair_interaction_nonbonded(it->first, it->second);
		}
		else {
			energy = _config_info->interaction->pair_interaction_term(_energy_term_id, it->first, it->second);
		}
		if(energy < _bond_threshold) {
			_bonds[it->first->index][it->second->index]++;
			_bonds[it->second->index][it->first->index]++;
		}
	}

	for(int i = 0; i < _config_info->N(); i++) {
		if(_bonds[i].size() > 0) {
			outstr << i << " ";
			for(map<int, int>::iterator it = _bonds[i].begin(); it != _bonds[i].end(); it++) {
				outstr << it->first << " ";
			}
			outstr << endl;
		}
	}

	if(interaction != nullptr) {
		interaction->no_three_body = false;
	}

	return outstr.str();
}
