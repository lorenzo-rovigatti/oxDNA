/*
 * PairEnergy.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: petr
 */

#include "PairEnergy.h"
//#include "../Interactions/RNAInteraction_V1.h"
#include <sstream>
#include <map>

template<typename number>
PairEnergy<number>::PairEnergy() {

}

template<typename number>
PairEnergy<number>::~PairEnergy() {

}

template<typename number>
void PairEnergy<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	int tmp = 0;
	tmp = 0;
	_print_all_particles = true;
	if (getInputInt(&my_inp,"particle1_id", &tmp, 0) == KEY_FOUND){
		_print_all_particles = false;
		_particle1_id = tmp;
		tmp = 0;
		getInputInt(&my_inp,"particle2_id", &tmp, 1);
		_particle2_id = tmp;
	}
}

template<typename number>
std::string PairEnergy<number>::get_output_string(llint curr_step) {
	BaseParticle<number> *p;
	BaseParticle<number> *q;


	std::stringstream output_str;

	number total_energy = 0.;
	number total_energy_diff = 0.;
	std::map<int, number> split_energies = this->_config_info.interaction->get_system_energy_split(this->_config_info.particles, *this->_config_info.N, this->_config_info.lists);

	if ((int)split_energies.size() == 7) output_str << "#id1 id2 FENE BEXC STCK NEXC HB CRSTCK CXSTCK total, t = " << curr_step << "\n";
	else output_str << "#id1 id2 FENE BEXC STCK NEXC HB CRSTCK CXSTCK DH total, t = " << curr_step << "\n";

	for(typename std::map<int,number>::iterator l = split_energies.begin(); l != split_energies.end(); l++)
		total_energy_diff += (*l).second;

	if (_print_all_particles){
		std::vector<ParticlePair<number> > neighbour_pairs = this->_config_info.lists->get_potential_interactions();

		//printf("Total fene energy is %f \n", split_energies[0]);
		for (int i = 0; i < (int)neighbour_pairs.size(); i++)
		{
			p = neighbour_pairs[i].first;
			q = neighbour_pairs[i].second;

			number pair_interaction = 0.;
			bool interaction_exists = false;
			std::stringstream pair_string;
			pair_string << p->index << " " << q->index;

			for(int k = 0; k < (int)split_energies.size(); k++)
			{
			   number k_th_interaction = this->_config_info.interaction->pair_interaction_term(k,q,p);
			   pair_interaction += k_th_interaction;
			   total_energy += k_th_interaction;
			   pair_string << " " <<  k_th_interaction ;
			   if(k_th_interaction != 0.f) interaction_exists = true;
			}
			pair_string << " " << pair_interaction << '\n';
			if (interaction_exists) output_str << pair_string.str();
		}
		output_str << "#Total energy per particle is  " << total_energy/ *this->_config_info.N << " and should be " << total_energy_diff / *this->_config_info.N << '\n';
	}
	else{
		p = this->_config_info.particles[_particle1_id];
		q = this->_config_info.particles[_particle2_id];

		number pair_interaction = 0.;
		std::stringstream pair_string;
		pair_string << p->index << " " << q->index;

		for(int k = 0; k < (int)split_energies.size(); k++)
		{
		   number k_th_interaction = this->_config_info.interaction->pair_interaction_term(k, q, p);
		   pair_interaction += k_th_interaction;
		   total_energy += k_th_interaction;
		   pair_string << " " <<  k_th_interaction ;
		}
		pair_string << " " << pair_interaction << '\n';
		// always print even if the interaction energy is zero
		output_str << pair_string.str();
		
	}

    return output_str.str();
}

template class PairEnergy<float>;
template class PairEnergy<double>;

