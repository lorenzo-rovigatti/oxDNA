/*
 * HBEnergy.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#include "HBEnergy.h"
#include "../Interactions/DNAInteraction.h"

HBEnergy::HBEnergy() {
	_mode = ALL_BASES;
}

HBEnergy::~HBEnergy() {

}

void HBEnergy::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	if(getInputString(&my_inp, "pairs_file", _list_file, 0) == KEY_FOUND) _mode = PAIRS_FROM_OP_FILE;
	if(getInputString(&my_inp, "bases_file", _list_file, 0) == KEY_FOUND) _mode = BASES_FROM_FILE;
}

void HBEnergy::init() {
	BaseObservable::init();

	switch(_mode) {
	case PAIRS_FROM_OP_FILE:
		_op.init_from_file(_list_file, _config_info->particles(), _config_info->N());
		break;
	case BASES_FROM_FILE: {
		ifstream inp(_list_file);
		while(!inp.eof()) {
			int n;
			inp >> n;
			_list.insert(n);
		}
		inp.close();
		break;
	}
	default:
		break;
	}
}

std::string HBEnergy::get_output_string(llint curr_step) {
	number energy = (number) 0.f;

	switch(_mode) {
	case PAIRS_FROM_OP_FILE: {
		vector_of_pairs inds = _op.get_hb_particle_list();
		for(vector_of_pairs::iterator i = inds.begin(); i != inds.end(); i++) {
			int p_ind = (*i).first;
			int q_ind = (*i).second;
			BaseParticle *p = _config_info->particles()[p_ind];
			BaseParticle *q = _config_info->particles()[q_ind];
			energy += _config_info->interaction->pair_interaction_term(DNAInteraction::HYDROGEN_BONDING, p, q);
		}
		energy *= (number) 2.f;
		break;
	}
	case BASES_FROM_FILE: {
		for(std::set<int>::iterator it = _list.begin(); it != _list.end(); it++) {
			BaseParticle *p = _config_info->particles()[*it];
			std::vector<BaseParticle *> neighs = _config_info->lists->get_all_neighbours(p);

			for(unsigned int j = 0; j < neighs.size(); j++) {
				energy += _config_info->interaction->pair_interaction_term(DNAInteraction::HYDROGEN_BONDING, p, neighs[j]);
			}
		}
		break;
	}
	default:
		energy = _config_info->interaction->get_system_energy_term(DNAInteraction::HYDROGEN_BONDING, _config_info->particles(), _config_info->lists);
		break;
	}

	return Utils::sformat("% 10.6lf", energy);
}
