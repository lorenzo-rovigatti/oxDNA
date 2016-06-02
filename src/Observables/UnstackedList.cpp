/*
 * StackeList.cpp
 *
 *  Created on: May 16, 2016
 *      Author: Ferdinando Randisi
 */

#include "UnstackedList.h"

#include <sstream>

template<typename number>
UnstackedList<number>::UnstackedList() {
	_threshold_fraction = 0.1;
	model = new Model;
}

template<typename number>
UnstackedList<number>::~UnstackedList() {
	delete model;
}

template<typename number> void
UnstackedList<number>::init(ConfigInfo<number> &config_info) {
   BaseObservable<number>::init(config_info);
}

template<typename number> void
UnstackedList<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	// Read temperature from the input file
  char TT[256];
  getInputString(&sim_inp, "T", TT, 1);    
  double T = Utils::get_temperature<number>(TT);
	// get the threshold fraction of the stacking energy
	getInputNumber(&my_inp,"threshold_fraction",&_threshold_fraction,0);
	// Make sure we're not using the sequence-dependent model
	bool average;
	if (getInputBool(&sim_inp,"use_average_seq",&average,0) == KEY_FOUND){
		if(!average){
		throw oxDNAException("%s: observable unstacked_list not implemented for the seqeunce dependent model",__FILE__);
		}
	}
	
	// Set the energy threshold
  getInputString(&sim_inp,"interaction_type",_interaction_type,1);
	number energy;
	if ( _interaction_type == "DNA"){
		energy = -(STCK_BASE_EPS_OXDNA + STCK_FACT_EPS_OXDNA * T);
	}
	else if ( _interaction_type == "DNA2" ){
		energy = -(STCK_BASE_EPS_OXDNA2 + STCK_FACT_EPS_OXDNA2 * T);
		
	}
	else if ( _interaction_type == "RNA" ){
		energy = -(model->RNA_STCK_BASE_EPS + model->RNA_STCK_FACT_EPS * T);
		// To implement, something like what has been done above has to be done
	}
	else throw oxDNAException("Interaction type %s not supported in observable unstacked_list (%s)",_interaction_type.c_str(),__FILE__);
		_threshold_energy = _threshold_fraction * energy;
	printf("%s interaction detected, threshold_energy is %f",_interaction_type.c_str(),_threshold_energy);
	
}

template<typename number>
std::string UnstackedList<number>::get_output_string(llint curr_step) {
	std::stringstream outstr;
	
		int N = *this->_config_info.N;
		for (int i=0;i<N-1;i++){
			BaseParticle<number> * p = this->_config_info.particles[i];
			BaseParticle<number> * q = this->_config_info.particles[i+1];

			number stacking_energy;
			if (_interaction_type == "DNA" || _interaction_type == "DNA2")
				stacking_energy = this->_config_info.interaction->pair_interaction_term(DNAInteraction<number>::STACKING,p,q);
			else if (_interaction_type == "RNA")
				stacking_energy = this->_config_info.interaction->pair_interaction_term(DNAInteraction<number>::STACKING,p,q);
			else throw oxDNAException("Interaction type %s not supported in observable unstacked_list (%s)",_interaction_type.c_str(),__FILE__);
			if (stacking_energy > _threshold_energy)
				outstr << i << " ";
		}
	return outstr.str();
	//TODO: debug & make sure it works, then keep carrying on with the plan
	/*
	// checking all particle pairs
	else{
		std::vector<ParticlePair<number> > pairs = this->_config_info.lists->get_potential_interactions();
		typename std::vector<ParticlePair<number> >::iterator it;
		for (it = pairs.begin(); it != pairs.end(); it ++ ) {
			BaseParticle<number> * p = (*it).first;
			BaseParticle<number> * q = (*it).second;
			int p_ind = p->get_index();
			int q_ind = q->get_index();
			number hb_energy = this->_config_info.interaction->pair_interaction_term(DNAInteraction<number>::HYDROGEN_BONDING, p, q);
			if (hb_energy < HB_CUTOFF){
				outstr << p_ind << " " << q_ind << "\n";
			}
		}
	}
	return outstr.str();
	*/
}

template class UnstackedList<float>;
template class UnstackedList<double>;
