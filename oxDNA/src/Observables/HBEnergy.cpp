/*
 * HBEnergy.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#include "HBEnergy.h"
#include "../Interactions/DNAInteraction.h"

template<typename number>
HBEnergy<number>::HBEnergy() {
	_mode = ALL_BASES;
}

template<typename number>
HBEnergy<number>::~HBEnergy() {

}

template<typename number>
void HBEnergy<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	if(getInputString(&my_inp, "pairs_file", _list_file, 0) == KEY_FOUND) _mode = PAIRS_FROM_OP_FILE;
	if(getInputString(&my_inp, "bases_file", _list_file, 0) == KEY_FOUND) _mode = BASES_FROM_FILE;
}

template<typename number>
void HBEnergy<number>::init(ConfigInfo<number> &config_info) {
   BaseObservable<number>::init(config_info);

   switch(_mode) {
   case PAIRS_FROM_OP_FILE:
	   _op.init_from_file(_list_file, this->_config_info.particles, *(this->_config_info.N));
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

template<typename number>
std::string HBEnergy<number>::get_output_string(llint curr_step) {
	number energy = (number) 0.f;

	switch(_mode) {
	case PAIRS_FROM_OP_FILE: {
		vector_of_pairs inds = _op.get_hb_particle_list();
		for (vector_of_pairs::iterator i = inds.begin(); i != inds.end(); i++) {
			int p_ind = (*i).first;
			int q_ind = (*i).second;
			BaseParticle<number> *p = this->_config_info.particles[p_ind];
			BaseParticle<number> *q = this->_config_info.particles[q_ind];
			energy += this->_config_info.interaction->pair_interaction_term(DNAInteraction<number>::HYDROGEN_BONDING, p, q);
		}
		energy *= (number) 2.f;
		break;
	}
	case BASES_FROM_FILE: {
		for(std::set<int>::iterator it = _list.begin(); it != _list.end(); it++) {
			BaseParticle<number> *p = this->_config_info.particles[*it];
			std::vector<BaseParticle<number> *> neighs = this->_config_info.lists->get_all_neighbours(p);

			for(unsigned int j = 0; j < neighs.size(); j++) {
				energy += this->_config_info.interaction->pair_interaction_term(DNAInteraction<number>::HYDROGEN_BONDING, p, neighs[j]);
			}
		}
		break;
	}
	default:
		energy = this->_config_info.interaction->get_system_energy_term(DNAInteraction<number>::HYDROGEN_BONDING, this->_config_info.particles, *this->_config_info.N, this->_config_info.lists);
		break;
	}

	return Utils::sformat("% 10.6lf", energy);
}

template class HBEnergy<float>;
template class HBEnergy<double>;
