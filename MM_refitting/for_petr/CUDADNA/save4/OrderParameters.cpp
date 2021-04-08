/*
 * OrderParameters.cpp
 *
 *  Created on: Feb 10, 2012
 *      Author: sulc
 */

#include "OrderParameters.h"
#include <cstdlib>

using namespace std;

HBParameter::HBParameter() {
	current_value = 0;
	stored_value = 0;
	name = "unnamed_parameter";
}

int HBParameter::plus_pair(base_pair& bp) {
	set_of_pairs::iterator i;
	i = counted_pairs.find(bp);
	if (i != counted_pairs.end()) {
		current_value++;
		return 1;
	} else {
		return 0;
	}
}

void HBParameter::save_pair_as_parameter(int a, int b) {
	base_pair bp(a, b);
	if (a > b) {
		bp.second = a;
		bp.first = b;
	}
	counted_pairs.insert(bp);
}

int HBParameter::minus_pair(base_pair& bp) {
	set_of_pairs::iterator i;
	i = counted_pairs.find(bp);
	if (i != counted_pairs.end()) {
		current_value--;
		if (current_value < 0) { //this should not happen!
			return -10;
		}
		return -1;
	} else {
		return 0;
	}
}

MinDistanceParameter::MinDistanceParameter() {
	current_value = -1;
	stored_value = -1;
	name = "unnamed_parameter";
}

void MinDistanceParameter::save_pair_as_parameter(int a, int b) {
	base_pair bp(a, b);
	if (a > b) {
		bp.second = a;
		bp.first = b;
	}
	counted_pairs.push_back(bp);
}

OrderParameters::OrderParameters() {
	_hb_parameters_count = 0;
	_distance_parameters_count = 0;
	_hb_states = NULL;
	_distance_states = NULL;
}

//adds given bonded pair to all values of order parameters
void OrderParameters::add_hb(int a, int b) {
	base_pair npair(a, b);
	if (a > b) {
		npair.first = b;
		npair.second = a;
	}
	for (int i = 0; i < _hb_parameters_count; i++) {
		_hb_parameters[i].plus_pair(npair);
	}
}

void OrderParameters::remove_hb(int a, int b) {
	base_pair npair(a, b);
	if (a > b) {
		npair.first = b;
		npair.second = a;
	}
	for (int i = 0; i < _hb_parameters_count; i++) {
		_hb_parameters[i].minus_pair(npair);
	}
}

void OrderParameters::reset(void) {
	for (int i = 0; i < _hb_parameters_count; i++) {
		_hb_parameters[i].reset();
	}
}

void OrderParameters::store(void) {
	for (int i = 0; i < _hb_parameters_count; i++) {
		_hb_parameters[i].store_current_value();
	}
	for (int i = 0; i < _distance_parameters_count; i++) {
		_distance_parameters[i].store_current_value();
	}
}

void OrderParameters::restore(void) {
	for (int i = 0; i < _hb_parameters_count; i++) {
		_hb_parameters[i].restore_current_value();
	}
	for (int i = 0; i < _distance_parameters_count; i++) {
		_distance_parameters[i].restore_current_value();
	}

}

int OrderParameters::get_hbpar_id_from_name(const char *name) {
	string stname(name);
	if (_hb_parnames.count(stname) != 0)
		return _hb_parnames[stname];
	else
		return -1;
}

int OrderParameters::get_distpar_id_from_name(const char *name) {
	string stname(name);

	if (_dist_parnames.count(stname) != 0)
		return _dist_parnames[stname];
	else
		return -1;
}

int *OrderParameters::get_hb_states(void) {
	for (int i = 0; i < _hb_parameters_count; i++)
		_hb_states[i] = _hb_parameters[i].get_state();

	return _hb_states;
}

int *OrderParameters::get_max_hb_states(void) {
	for (int i = 0; i < _hb_parameters_count; i++)
		_hb_states[i] = _hb_parameters[i].get_max_state();

	return _hb_states;
}

double * OrderParameters::get_distance_states(void) {
	for (int i = 0; i < _distance_parameters_count; i++)
		_distance_states[i] = _distance_parameters[i].get_state();

	return _distance_states;
}

OrderParameters::~OrderParameters() {
	delete[] _hb_states;
	delete[] _distance_states;
}

