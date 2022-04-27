/*
 * OrderParameterValues.cpp
 *
 *  Created on: Mar 2, 2013
 *      Author: petr
 */

#include "OrderParameterValues.h"

#include <sstream>

OrderParameterValues::OrderParameterValues() {

}

OrderParameterValues::~OrderParameterValues() {

}

void OrderParameterValues::init() {
	BaseObservable::init();

	ifstream tmpf(_order_parameters_file);
	if(!tmpf.good()) throw oxDNAException("(OrderParameterValues.cpp) Can't read file '%s'", _order_parameters_file);

	_op.init_from_file(_order_parameters_file, _config_info->particles(), _config_info->N());

}

void OrderParameterValues::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	bool parameters_loaded = false;
	//first try load parameters from specific op_file key
	if(getInputString(&my_inp, "order_parameters_file", _order_parameters_file, 0) == KEY_NOT_FOUND) {
		if(getInputString(&my_inp, "op_file", _order_parameters_file, 0) != KEY_NOT_FOUND) {
			parameters_loaded = true;
		}

	}
	else {
		parameters_loaded = true;
	}

	//no success, so try it load from the main simulation input file
	if(parameters_loaded == false) {
		if(getInputString(&sim_inp, "order_parameters_file", _order_parameters_file, 0) == KEY_NOT_FOUND) {
			getInputString(&sim_inp, "op_file", _order_parameters_file, 1);
		}
	}

}

std::string OrderParameterValues::get_output_string(llint curr_step) {
	_op.reset();
	_op.fill_distance_parameters(_config_info->particles(), _config_info->box);

	std::vector<ParticlePair> neighbour_pairs = _config_info->lists->get_potential_interactions();

	for(int i = 0; i < (int) neighbour_pairs.size(); i++) {
		BaseParticle * p = neighbour_pairs[i].first;
		BaseParticle * q = neighbour_pairs[i].second;

		number hb_energy = _config_info->interaction->pair_interaction_term(DNAInteraction::HYDROGEN_BONDING, p, q);
		if(hb_energy < HB_CUTOFF) _op.add_hb(p->index, q->index);
	}

	std::stringstream outstr;
	for(int i = 0; i < _op.get_hb_parameters_count(); i++) {
		outstr << _op.get_hb_parameter(i) << " ";
	}
	for(int i = 0; i < _op.get_distance_parameters_count(); i++) {
		outstr << _op.get_distance_parameter(i) << " ";
	}

	return outstr.str();
}
