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
std::string PairEnergy<number>::get_output_string(llint curr_step) {
	BaseParticle<number> *p;
	BaseParticle<number> *q;


	std::stringstream output_str;
	output_str << "#id1 id2 FENE BEXC STCK NEXC HB CRSTCK CXSTCK total, t = " << curr_step << "\n";

	number total_energy = 0.;
	number total_energy_diff = 0.;
	//std::map<int,number> split_energies = this->_config_info->get_system_energy_split(this->_config_info.particles, *this->_config_info.N, *this->_config_info.box_side);
	std::map<int, number> split_energies = this->_config_info.interaction->get_system_energy_split(this->_config_info.particles, *this->_config_info.N, *this->_config_info.box_side);

	for(typename std::map<int,number>::iterator l = split_energies.begin(); l != split_energies.end(); l++)
		total_energy_diff += (*l).second;

	std::vector<std::pair<BaseParticle<number> *, BaseParticle<number> *> > neighbour_pairs = this->_config_info.interaction->get_potential_interactions(this->_config_info.particles, *this->_config_info.N, *this->_config_info.box_side);

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

	//number energy = this->_config_info.interaction->get_system_energy_term(RNAInteraction<number>::HYDROGEN_BONDING, this->_config_info.particles, *this->_config_info.N, *this->_config_info.box_side);
	//energy /= *this->_config_info.N;
    output_str << "#Total energy per particle is  " << total_energy/ *this->_config_info.N << " and should be " << total_energy_diff / *this->_config_info.N << '\n';
    //printf("Finished, I have %s",output_str.str().c_str());
    return output_str.str();
}

template class PairEnergy<float>;
template class PairEnergy<double>;

