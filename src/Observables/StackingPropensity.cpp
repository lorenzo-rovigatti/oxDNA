/*
 * StackingPropensity.cpp
 *
 *  Created on: Feb 13, 2023
 *      Author: andrea_bonato
 */

#include "StackingPropensity.h"
#include "../Utilities/Utils.h"
#include "../Utilities/oxDNAException.h"
#include "../model.h"
#include <sstream>
#include <map>
#include <math.h>

StackingPropensity::StackingPropensity() {
	_print_header = true;
	_print_all_particles = false;
	_bonded_only = true;
	_average = true;
	_first = true;
}

StackingPropensity::~StackingPropensity() {

}

void StackingPropensity::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputBool(&my_inp, "print_all", &_print_all_particles, 0);
	getInputBool(&my_inp, "print_header", &_print_header, 0);
	getInputBool(&my_inp, "bonded_only", &_bonded_only, 0);
	
	// Read interaction type
	std::string inter_type("");
	getInputString(&sim_inp, "interaction_type", inter_type, 0);
	
	oxdna_version = -1;

	if(inter_type.compare("DNA") == 0) {
		 oxdna_version = 1;
	}
	else if(inter_type.compare("DNA2") == 0) {
		 oxdna_version = 2;
	}
	else if(inter_type.compare("DNA2SD") == 0) {
		 oxdna_version = 3;
	}
	else if(inter_type.compare("DNA3") == 0) {
		 oxdna_version = 3;
	}
	
	char T[256];
	getInputString(&sim_inp, "T", T, 1);
	number _T = Utils::get_temperature(T); //get temperature
	
	// set the default values
	for(int i = 0; i < 5; i++) {
		for(int j = 0; j < 5; j++) {
			// stacking
			if(oxdna_version == 1) min_stacking_energy[i][j] = -1.0*(STCK_BASE_EPS_OXDNA + STCK_FACT_EPS_OXDNA * _T)*SQR(1 - exp(-(STCK_RC - STCK_R0) * STCK_A)); //XXXTODO
			else if(oxdna_version == 2) min_stacking_energy[i][j] = -1.0 * (STCK_BASE_EPS_OXDNA2 + STCK_FACT_EPS_OXDNA2 * _T) * SQR(1 - exp(-(STCK_RC - STCK_R0) * STCK_A));
			else if(oxdna_version == 3) min_stacking_energy[i][j] = -1.0 * (STCK_BASE_EPS_OXDNA2 + STCK_FACT_EPS_OXDNA2 * _T) * SQR(1 - exp(-(STCK_RC - STCK_R0) * STCK_A)); //XXXTODO
			else min_stacking_energy[i][j] = 0.;
		}
	}
	
	//read sequence specific file (if sequence specific interactions enabled)
	if(getInputBool(&sim_inp, "use_average_seq", &_average, 0) == KEY_FOUND) {
		if(!_average) {
			getInputString(&sim_inp, "seq_dep_file", _seq_filename, 1);
			OX_LOG(Logger::LOG_INFO, "Using '%s' as the input for sequence-dependent values", _seq_filename.c_str());
		}
	}
		
	// if we want to use the sequence dependence then we overwrite all the
	// epsilon regarding interactions between true bases
	if(!_average) {
		char key[256];
		float tmp_value, stck_fact_eps;

		input_file seq_file;
		seq_file.init_from_filename(_seq_filename.c_str());
		if(seq_file.state == ERROR)
			throw oxDNAException("Caught an error while opening sequence dependence file '%s'", _seq_filename.c_str());

		// stacking
		
		getInputFloat(&seq_file, "STCK_FACT_EPS", &stck_fact_eps, 1);
		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < 4; j++) {
				sprintf(key, "STCK_%c_%c", Utils::encode_base(i), Utils::encode_base(j));
				getInputFloat(&seq_file, key, &tmp_value, 1);
				min_stacking_energy[i][j] = -1.0 * tmp_value * (1.0 - stck_fact_eps + (_T * 9.0 * stck_fact_eps)) * SQR(1 - exp(-(STCK_RC - STCK_R0) * STCK_A));
			}
		}
		
	}
}

std::string StackingPropensity::get_output_string(llint curr_step) {

	if (oxdna_version == -1) return "";

	std::stringstream output_str;

	//number total_energy = 0.;
	//number total_energy_diff = 0.;
	//std::map<int, number> split_energies = _config_info->interaction->get_system_energy_split(_config_info->particles(), _config_info->lists);


	
	if(_first && _print_header) {
		
		if(!_print_all_particles) output_str << "#time Npairs Nstacked\n";
		else output_str << "#bp_id base1_id base2_id  stk_energy stacked(1=yes; energy > 0.3 min) \n";
		//else output_str << "#t = " << curr_step << "\n";
			if(!_average) {
		output_str << "#min_stk_energies stacked_threshold: \n";
		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < 4; j++) {
				output_str << Utils::encode_base(i) << Utils::encode_base(j) << ": " <<  min_stacking_energy[i][j] << "  " << 0.3*min_stacking_energy[i][j] << "\n";
				}
			}
			output_str << "\n";	
		}
		else {
			output_str << "#min_stk_energy: " << min_stacking_energy[0][0] << "\n\n";
		}
		
		_first = false;
	}
	

	if(_print_all_particles) output_str << curr_step << "\n";
	
	std::vector<ParticlePair> neigh_pairs; 

	if(_bonded_only) {
		for(auto p : _config_info->particles()) {
			for(unsigned int n = 0; n < p->affected.size(); n++) {
				BaseParticle *p1 = p->affected[n].first;
				BaseParticle *p2 = p->affected[n].second;
				if(p1->index <= p->index && p2->index <= p->index) {
					neigh_pairs.push_back(ParticlePair(p1, p2));
				}
			}
		}
	}
	
	
	else {
	
		neigh_pairs = _config_info->lists->get_potential_interactions(); //list of interacting pairs (bonded and not bonded)
	}


	BaseParticle *p;
	BaseParticle *q;
	
	int STCK_N = 0;
	int N = 0;
	for(unsigned int i = 0; i < neigh_pairs.size(); i++) {
		p = neigh_pairs[i].first;
		q = neigh_pairs[i].second;

		//number pair_interaction = 0.;
		//bool interaction_exists = false;
		std::stringstream pair_string;
		pair_string << p->index << " " << q->index;

		int k = 2;  //STACKING, see DNAInteracion.h
		number stacking_energy = _config_info->interaction->pair_interaction_term(k, q, p);
		pair_string << " " << stacking_energy;
		
		if(stacking_energy < 0.3*min_stacking_energy[p->type][q->type]) {
			STCK_N++;
			pair_string << " 1\n";
		}
		else pair_string << " 0\n";		
				
		if(_print_all_particles) output_str << pair_string.str();
		N++;
	}
	if(!_print_all_particles) output_str << curr_step << " " << N << " " << STCK_N;
	
	return output_str.str();
	
}

