/*
 * PLCluster.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: petr
 */

#include "PLCluster.h"
#include "../Interactions/PatchyShapeInteraction.h"
//#include "../Particles/PatchyShapeParticle.h"

#include <sstream>
#include <map>

template<typename number>
PLCluster<number>::PLCluster() {
 _show_types = true;
}

template<typename number>
PLCluster<number>::~PLCluster() {

}

template<typename number>
void PLCluster<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	int show_types = 0;

	if( getInputBoolAsInt(&my_inp,"show_types",&show_types,1) == KEY_FOUND)
	{
		this->_show_types = (bool)show_types;
	}

}

template<typename number>
std::string PLCluster<number>::get_output_string(llint curr_step) {

	PatchyShapeInteraction<number> *interaction = dynamic_cast<PatchyShapeInteraction<number> * >(this->_config_info.interaction);
    interaction->check_patchy_locks();

	BaseParticle<number> *p;
	int N = *this->_config_info.N;
	int *cluster_membership = new int[N];
	for (int i = 0; i < N; i++)
		cluster_membership[i] = -1;
	//BaseParticle<number> *q;
	number pair_energy;
	//number total_energy;
	int cluster_index = 0;
	std::stringstream output_str;
	int cluster_count = 0;

	for (int i = 0; i < *this->_config_info.N; i++) {
		p = this->_config_info.particles[i];
		std::vector<BaseParticle<number> *> neighs =
				this->_config_info.lists->get_all_neighbours(p);

        //printf("Particle %d has %d neighbors \n",i,neighs.size());
		for (unsigned int j = 0; j < neighs.size(); j++) {
			pair_energy = this->_config_info.interaction->pair_interaction_term(
					PatchyShapeInteraction<number>::PATCHY, p, neighs[j]);
			if (pair_energy < 0 && i < neighs[j]->index) {
				if (cluster_membership[p->index] == -1
						&& cluster_membership[neighs[j]->index] == -1) {
					cluster_membership[p->index] =
							cluster_membership[neighs[j]->index] =
									cluster_index;
					cluster_index++;
				} else if (cluster_membership[p->index] == -1
						&& cluster_membership[neighs[j]->index] > -1) {
					cluster_membership[p->index] =
							cluster_membership[neighs[j]->index];
				} else if (cluster_membership[p->index] > -1
						&& cluster_membership[neighs[j]->index] == -1) {
					cluster_membership[neighs[j]->index] =
							cluster_membership[p->index];
				} else if (cluster_membership[p->index] > -1
						&& cluster_membership[neighs[j]->index] > -1
						&& cluster_membership[neighs[j]->index]
								!= cluster_membership[p->index]) {
					int oldindex = cluster_membership[neighs[j]->index];
					cluster_membership[neighs[j]->index] =
							cluster_membership[p->index];
					for (int k = 0; k < N; k++) {
						if (cluster_membership[k] == oldindex) {
							cluster_membership[k] =
									cluster_membership[p->index];
						}
					}
				}
			}
		}


	}

	if (cluster_index > 0) {
		    vector<int> cluster_members;
			for (int j = 0; j < cluster_index; j++) {
				cluster_members.clear();
				std::stringstream clusout;
				clusout << "( ";
				bool isempty = true;
				for (int k = 0; k < N; k++) {
					if (cluster_membership[k] == j) {
						if(! this->_show_types)
						{
						     //clusout << k << " ";
						     cluster_members.push_back(k);
						}
						else
						{
							//clusout << this->_config_info.particles[k]->type << " " ;
							cluster_members.push_back( this->_config_info.particles[k]->type );
						}
						isempty = false;
					}
				}
				//clusout << ")";
				if (!isempty) {
					cluster_count++;
					std::sort(cluster_members.begin(),cluster_members.end());
				    for(vector<int>::iterator k = cluster_members.begin(); k != cluster_members.end(); ++k)
					{
					  clusout << *k << " ";
					}
				    clusout << ")";
					output_str << clusout.str() << " ";
				}
			}
	}


	delete [] cluster_membership;
	std::stringstream out;
	out << cluster_count << " " << output_str.str();
	return out.str();
}

template class PLCluster<float>;
template class PLCluster<double>;

