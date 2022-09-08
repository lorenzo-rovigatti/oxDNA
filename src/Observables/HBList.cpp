/*
 * HBList.cpp
 *
 *  Created on: Apr 11, 2013
 *      Author: Ben Snodin, later modified by Ferdinando in e.g. September 2017.
 */

#include "HBList.h"

#include <sstream>

HBList::HBList() {
	_max_shift = 10;
	_measure_mean_shift = false;
	_only_count = false;
	_read_op = true;
	_mean_shift = 0.;
}

HBList::~HBList() {

}

void HBList::init() {
	BaseObservable::init();
	if(_read_op) {
		_op.init_from_file(_order_parameters_file, _config_info->particles(), _config_info->N());
	}
}

void HBList::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	bool parameters_loaded = false;
	//first try to load parameters from specific op_file key
	if(getInputString(&my_inp, "order_parameters_file", _order_parameters_file, 0) == KEY_NOT_FOUND) {
		if(getInputString(&my_inp, "op_file", _order_parameters_file, 0) != KEY_NOT_FOUND) {
			parameters_loaded = true;
		}

	}
	else {
		parameters_loaded = true;
	}

	//no success, so try to load it from the main simulation input file
	if(parameters_loaded == false) {
		if(getInputString(&sim_inp, "order_parameters_file", _order_parameters_file, 0) == KEY_NOT_FOUND) {
			if(getInputString(&sim_inp, "op_file", _order_parameters_file, 0) == KEY_NOT_FOUND) {
				OX_LOG(Logger::LOG_INFO, "No order parameters file specified in input file; printing indices of any particle pairs that have a hydrogen bond between them");
				_read_op = false;
			}
		}
	}
	// read the _only_count and _measure_mean_shift flags
	getInputBool(&my_inp, "only_count", &_only_count, 0);
	if(getInputBool(&my_inp, "measure_mean_shift", &_measure_mean_shift, 0) == KEY_FOUND) {
		getInputInt(&my_inp, "max_shift", &_max_shift, 0);
	}

}

std::vector<std::pair<int, int>> HBList::hb_list() {
	std::vector<std::pair<int, int>> result;

	unsigned long long int total_reg = 0;
	unsigned long long int n_reg_entries = 0;

	// using an order parameters file
	if(_read_op) {
		vector_of_pairs inds = _op.get_hb_particle_list();
		for(vector_of_pairs::iterator i = inds.begin(); i != inds.end(); i++) {
			int p_ind = (*i).first;
			int q_ind = (*i).second;
			BaseParticle *p = _config_info->particles()[p_ind];
			BaseParticle *q = _config_info->particles()[q_ind];
			number hb_energy = _config_info->interaction->pair_interaction_term(DNAInteraction::HYDROGEN_BONDING, p, q);

			if(hb_energy < HB_CUTOFF) {
				result.emplace_back(p_ind, q_ind);
			}
		}
	}
	// checking all particle pairs
	else {
		std::vector<ParticlePair> pairs = _config_info->lists->get_potential_interactions();
		typename std::vector<ParticlePair>::iterator it;
		for(it = pairs.begin(); it != pairs.end(); it++) {
			BaseParticle * p = (*it).first;
			BaseParticle * q = (*it).second;
			int p_ind = p->get_index();
			int q_ind = q->get_index();
			number hb_energy = _config_info->interaction->pair_interaction_term(DNAInteraction::HYDROGEN_BONDING, p, q);
			// what to do if it's unbound
			int complementary_ind = _config_info->N() - 1 - p_ind;
			if(hb_energy < HB_CUTOFF) {
				result.emplace_back(p_ind, q_ind);
			}
			else if(_measure_mean_shift && complementary_ind == q_ind) {
				//this always takes into account periodic boundary conditions, which makes sense.
				LR_vector distance_vector = _config_info->box->min_image(p->pos, q->pos);
				double distance = distance_vector * distance_vector;
				BaseParticle * ahead = q;
				BaseParticle * before = q;
				//find the register for this particular base-pair
				int reg = 0;
				for(int i = 0; i < _max_shift; i++) {
					double newdist = -1;
					if(ahead->n5 != P_VIRTUAL) {
						ahead = ahead->n5;
						LR_vector newdist_vector = _config_info->box->min_image(p->pos, ahead->pos);
						newdist = newdist_vector * newdist_vector;
						if(newdist < distance) {
							distance = newdist;
							reg = i + 1;
						}
					}
					if(before->n3 != P_VIRTUAL) {
						before = before->n3;
						LR_vector newdist_vector = _config_info->box->min_image(p->pos, before->pos);
						newdist = newdist_vector * newdist_vector;
						if(newdist < distance) {
							distance = newdist;
							reg = i + 1;
						}
					}
				}
				total_reg += reg;
				n_reg_entries++;
			}
		}
	}

	_mean_shift = (n_reg_entries > 0) ? total_reg / (number) n_reg_entries : (number) -1.0;

	return result;
}

std::string HBList::get_output_string(llint curr_step) {
	std::stringstream outstr;
	auto list = hb_list();

	if(!_only_count and !_measure_mean_shift) {
		outstr << "# step " << curr_step << std::endl;
		for(auto &hb_pair : list) {
			outstr << hb_pair.first << " " << hb_pair.second << std::endl;
		}
	}
	if(_only_count) {
		outstr << list.size();
	}
	else if(_measure_mean_shift) {
		if(_mean_shift == (number) -1.0) {
			outstr << "none";
		}
		else {
			outstr << _mean_shift;
		}
	}

	return outstr.str();
}
