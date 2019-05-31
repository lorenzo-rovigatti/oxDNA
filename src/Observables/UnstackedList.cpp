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
	input_file seq_file;
	if(getInputBool(&sim_inp,"use_average_seq",&average,0) == KEY_FOUND){
		if(!average){
			std::string seq_filename;
			getInputString(&sim_inp, "seq_dep_file", seq_filename, 1);
			loadInputFile(&seq_file, seq_filename.c_str());
		}
	}

	// Set the energy threshold
  getInputString(&sim_inp,"interaction_type",_interaction_type,1);
	number energy;
	if(_interaction_type == "DNA" or _interaction_type == "DNA2" or _interaction_type == "DNA2ModInteraction" or _interaction_type == "DNA2_nomesh"){
		for(int i = 0; i < 5; i++) {
			for(int j = 0; j < 5; j++) {
				if(_interaction_type == "DNA" or _interaction_type == "DNA_nomesh")
					energy = -(STCK_BASE_EPS_OXDNA + STCK_FACT_EPS_OXDNA * T);
				else if(_interaction_type == "DNA2" or _interaction_type == "DNA2ModInteraction" or _interaction_type == "DNA2_nomesh")
					energy = -(STCK_BASE_EPS_OXDNA2 + STCK_FACT_EPS_OXDNA2 * T);
				else throw oxDNAException("In observable UnstackedList.cpp: don't know what to do with interaction \"%s\"",_interaction_type.c_str());
				_threshold_energies[i][j] = energy * _threshold_fraction;
			}
		}
		if(!average){
			char key[256];
			float tmp_value, stck_fact_eps;
			getInputFloat(&seq_file, "STCK_FACT_EPS", &stck_fact_eps, 1);
			for(int i = 0; i < 5; i++) {
				for(int j = 0; j < 5; j++) {
					sprintf(key, "STCK_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
					if(getInputFloat(&seq_file, key, &tmp_value, !(i == 4 or j == 4)) == KEY_FOUND){//interactions that involve dummy bases are not mandatory
						energy =  - tmp_value * ( 1.0 - stck_fact_eps + (T*9.0 * stck_fact_eps ) ); 
						_threshold_energies[i][j] = energy * _threshold_fraction;
					}
				}
			}
		}
	}
	else if(_interaction_type == "RNA" or _interaction_type == "RNA2"){
		for(int i = 0; i < 5; i++) {
			for(int j = 0; j < 5; j++) {
				// stacking
				energy = -(model->RNA_STCK_BASE_EPS + model->RNA_STCK_FACT_EPS * T);
				_threshold_energies[i][j] = energy * _threshold_fraction;
			}
		}
		if(!average){	
			// stacking
			char key[256];
			float tmp_value, stck_fact_eps;
			getInputFloat(&seq_file, "ST_T_DEP", &stck_fact_eps, 1);
			for(int i = 0; i < 4; i++) {
				for(int j = 0; j < 4; j++) {
					sprintf(key, "STCK_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
					getInputFloat(&seq_file, key, &tmp_value, 1);
					energy =  -tmp_value * ( 1.0 + T*stck_fact_eps  );
					_threshold_energies[i][j] = energy * _threshold_fraction;
				}
			}
		}
	}
	else throw oxDNAException("Interaction type %s not supported in observable unstacked_list (%s)",_interaction_type.c_str(),__FILE__);
}

template<typename number>
std::string UnstackedList<number>::get_output_string(llint curr_step){
	std::stringstream outstr;

	int N = *this->_config_info.N;
	for(int i=0;i<N-1;i++){
		BaseParticle<number> * p = this->_config_info.particles[i];
		BaseParticle<number> * q = this->_config_info.particles[i+1];

		//the following line don't need to be different for every interaction, assuming that all these interactions have the same values of STACKING and COAXIAL_STACKING, which is true
		// for any interaction that inherits from DNAInteraction (unless it's changed explicitly).
		number stacking_energy = this->_config_info.interaction->pair_interaction_term(DNAInteraction<number>::STACKING,p,q) + this->_config_info.interaction->pair_interaction_term(DNAInteraction<number>::COAXIAL_STACKING,p,q);
		if(stacking_energy > _threshold_energies[p->type][q->type])
			outstr << i << " ";
	}
	return outstr.str();
}

template class UnstackedList<float>;
template class UnstackedList<double>;
