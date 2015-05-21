/*
 * MCBackend.cpp
 *
 *  Created on: 25/nov/2010
 *      Author: lorenzo
 */

#include <sstream>

#include "MCBackend.h"

template<typename number>
MCBackend<number>::MCBackend() : SimBackend<number>(), _MC_moves(3), _overlap(false), _check_energy_counter(0) {
	this->_sim_type = SIM_MC;

	_delta = new number[_MC_moves];
	_tries = new llint[_MC_moves];
	_accepted = new llint[_MC_moves];
	_adjust_moves = false;
	_MC_equilibration_steps = 0;
	_check_energy_every = 10;

	for(int i = 0; i < _MC_moves; i++) _tries[i] = _accepted[i] = 0;
}

template<typename number>
MCBackend<number>::~MCBackend() {
	delete[] _delta;
	delete[] _tries;
	delete[] _accepted;
}

template<>
void MCBackend<float>::_get_number_settings(input_file &inp) {
	if(getInputFloat(&inp, "check_energy_threshold", &_check_energy_threshold, 0) == KEY_NOT_FOUND)
		_check_energy_threshold = 1e-2f;
	getInputFloat(&inp, "delta_translation", &_delta[MC_MOVE_TRANSLATION], 1);
	getInputFloat(&inp, "delta_rotation", &_delta[MC_MOVE_ROTATION], 1);
	if(_ensemble == MC_ENSEMBLE_NPT) getInputFloat(&inp, "delta_volume", &_delta[MC_MOVE_VOLUME], 1);
	if(_ensemble == MC_ENSEMBLE_NPT) getInputFloat(&inp, "P", &(this->_P), 1);
}

template<>
void MCBackend<double>::_get_number_settings(input_file &inp) {
	if(getInputDouble(&inp, "check_energy_threshold", &_check_energy_threshold, 0) == KEY_NOT_FOUND)
		_check_energy_threshold = 1e-6;
	getInputDouble(&inp, "delta_translation", &_delta[MC_MOVE_TRANSLATION], 1);
	getInputDouble(&inp, "delta_rotation", &_delta[MC_MOVE_ROTATION], 1);
	if(_ensemble == MC_ENSEMBLE_NPT) getInputDouble(&inp, "delta_volume", &_delta[MC_MOVE_VOLUME], 1);
	if(_ensemble == MC_ENSEMBLE_NPT) getInputDouble(&inp, "P", &(this->_P), 1);
}

template<typename number>
void MCBackend<number>::get_settings(input_file &inp) {
	char tmp[256];

	SimBackend<number>::get_settings(inp);

	getInputInt(&inp, "check_energy_every", &_check_energy_every, 0);
	getInputString(&inp, "ensemble", tmp, 1);
	if(strncasecmp(tmp, "npt", 256) == 0) {
		_ensemble = MC_ENSEMBLE_NPT;
		//throw oxDNAException("Ensemble NPT will be implemented soon");
		int check = getInputString(&inp, "list_type", tmp, 0);
		if (check == KEY_NOT_FOUND || (strncasecmp(tmp, "cells", 256) != 0 && strncasecmp(tmp, "no", 256) != 0)) throw oxDNAException ("NPT ensemble requires no lists or cells to handle interaction lists; set list_type=cells in the input file");

	}
	else if(strncasecmp(tmp, "nvt", 256) == 0) _ensemble = MC_ENSEMBLE_NVT;
	else throw oxDNAException("Ensemble '%s' not supported\n", tmp);

	int tmpi;
	if (getInputBoolAsInt (&inp, "adjust_moves", &tmpi, 0) == KEY_FOUND) {
		_adjust_moves = (tmpi > 0);
		if (_adjust_moves) OX_LOG(Logger::LOG_INFO, "(MCBackend) adjusting moves in the equilibration phase");
		getInputLLInt(&inp, "equilibration_steps", &_MC_equilibration_steps, 1);
	}

	_get_number_settings(inp);

	char energy_file[512];
	llint print_every;
	getInputString(&inp, "energy_file", energy_file, 1);
	getInputLLInt(&inp, "print_energy_every", &print_every, 1);
	// we build the default stream of observables;
	// we build a helper string for that
	std::string fake = Utils::sformat("{\n\tname = %s\n\tprint_every = %lld\n}\n", energy_file, print_every);
	this->_obs_output_file = new ObservableOutput<number>(fake, inp);
	this->_obs_output_file->add_observable("type = step");
	this->_obs_output_file->add_observable("type = potential_energy");
	if (_ensemble == MC_ENSEMBLE_NPT) this->_obs_output_file->add_observable("type = density");
	this->_obs_output_file->add_observable("type = backend_info");
	this->_obs_outputs.push_back(this->_obs_output_file);

	// now we do the same thing for stdout
	int no_stdout_energy = 0;
	getInputBoolAsInt(&inp, "no_stdout_energy", &no_stdout_energy, 0);
	if(!no_stdout_energy) {
		fake = Utils::sformat("{\n\tname = %s\n\tprint_every = %lld\n}\n", "stdout", print_every);
		this->_obs_output_stdout = new ObservableOutput<number>(fake, inp);
		this->_obs_outputs.push_back(this->_obs_output_stdout);
		this->_obs_output_stdout->add_observable("type = step");
		this->_obs_output_stdout->add_observable("type = potential_energy");
		if (_ensemble == MC_ENSEMBLE_NPT) this->_obs_output_stdout->add_observable("type = density");
		this->_obs_output_stdout->add_observable("type = backend_info");
	}
}

template<typename number>
void MCBackend<number>::print_observables(llint curr_step) {
	std::string tmpstr("");
	for(int i = 0; i < _MC_moves; i++) {
		number ratio = (this->_tries[i] > 0) ? this->_accepted[i]/(float)this->_tries[i] : 0;
		//this->_backend_info += Utils::sformat("  %5.3lf", ratio);
		tmpstr += Utils::sformat("  %5.3lf", ratio);
	}
	this->_backend_info.insert(0, tmpstr + "  ");

	SimBackend<number>::print_observables(curr_step);
}

template class MCBackend<float>;
template class MCBackend<double>;
