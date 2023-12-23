/*
 * MCBackend.cpp
 *
 *  Created on: 25/nov/2010
 *      Author: lorenzo
 */

#include "MCBackend.h"

#include "../Observables/ObservableOutput.h"
#include <sstream>

MCBackend::MCBackend() :
				SimBackend(),
				_MC_moves(3),
				_overlap(false),
				_check_energy_counter(0) {
	_delta.resize(_MC_moves);
	_tries.resize(_MC_moves);
	_accepted.resize(_MC_moves);
	_adjust_moves = false;
	_MC_equilibration_steps = 0;
	_check_energy_every = 10;
	_ensemble = -1;
	_check_energy_threshold = 1e-6;

	for(int i = 0; i < _MC_moves; i++)
		_tries[i] = _accepted[i] = 0;
}

MCBackend::~MCBackend() {

}

void MCBackend::get_settings(input_file &inp) {
	char tmp[256];

	SimBackend::get_settings(inp);

	std::string sim_type;
	getInputString(&inp, "sim_type", sim_type, 0);

	getInputInt(&inp, "check_energy_every", &_check_energy_every, 0);
	if(getInputNumber(&inp, "check_energy_threshold", &_check_energy_threshold, 0) == KEY_NOT_FOUND) {
		_check_energy_threshold = 1e-6;
	}

	if(sim_type == "MC") {
		getInputString(&inp, "ensemble", tmp, 1);
		if(strncasecmp(tmp, "npt", 256) == 0) {
			_ensemble = MC_ENSEMBLE_NPT;
			int check = getInputString(&inp, "list_type", tmp, 0);
			if(check == KEY_NOT_FOUND || (strncasecmp(tmp, "cells", 256) != 0 && strncasecmp(tmp, "no", 256) != 0 && strncasecmp(tmp, "rodcells", 256) != 0))
				throw oxDNAException("NPT ensemble requires no lists or cells to handle interaction lists; set list_type=cells in the input file");

		}
		else if(strncasecmp(tmp, "nvt", 256) == 0) {
			_ensemble = MC_ENSEMBLE_NVT;
		}
		else {
			throw oxDNAException("Ensemble '%s' not supported\n", tmp);
		}

		int tmpi;
		if(getInputBoolAsInt(&inp, "adjust_moves", &tmpi, 0) == KEY_FOUND) {
			_adjust_moves = (tmpi > 0);
			if(_adjust_moves)
				OX_LOG(Logger::LOG_INFO, "(MCBackend) adjusting moves in the equilibration phase");
			getInputLLInt(&inp, "equilibration_steps", &_MC_equilibration_steps, 1);
		}

		getInputNumber(&inp, "delta_translation", &_delta[MC_MOVE_TRANSLATION], 1);
		getInputNumber(&inp, "delta_rotation", &_delta[MC_MOVE_ROTATION], 1);

		if(_ensemble == MC_ENSEMBLE_NPT) {
			getInputNumber(&inp, "P", &(_P), 1);
			getInputNumber(&inp, "delta_volume", &_delta[MC_MOVE_VOLUME], 1);
		}
	}
	else if(sim_type == "VMMC" || sim_type == "PT_VMMC") {
		getInputNumber(&inp, "delta_translation", &_delta[MC_MOVE_TRANSLATION], 1);
		getInputNumber(&inp, "delta_rotation", &_delta[MC_MOVE_ROTATION], 1);
	}

	char energy_file[512];
	llint print_every;
	getInputString(&inp, "energy_file", energy_file, 1);
	getInputLLInt(&inp, "print_energy_every", &print_every, 1);
	// we build the default stream of observables;
	// we build a helper string for that
	std::string fake = Utils::sformat("{\n\tname = %s\n\tprint_every = %lld\n}\n", energy_file, print_every);
	_obs_output_file = std::make_shared<ObservableOutput>(fake);
	_obs_output_file->add_observable("type = step");
	_obs_output_file->add_observable("type = potential_energy");
	if(_ensemble == MC_ENSEMBLE_NPT) {
		_obs_output_file->add_observable("type = density");
	}
	_obs_output_file->add_observable("type = backend_info");
	add_output(_obs_output_file);

	// now we do the same thing for stdout
	int no_stdout_energy = 0;
	getInputBoolAsInt(&inp, "no_stdout_energy", &no_stdout_energy, 0);
	if(!no_stdout_energy) {
		fake = Utils::sformat("{\n\tname = %s\n\tprint_every = %lld\n}\n", "stdout", print_every);
		_obs_output_stdout = std::make_shared<ObservableOutput>(fake);
		add_output(_obs_output_stdout);
		_obs_output_stdout->add_observable("type = step");
		_obs_output_stdout->add_observable("type = potential_energy");
		if(_ensemble == MC_ENSEMBLE_NPT) {
			_obs_output_stdout->add_observable("type = density");
		}
		_obs_output_stdout->add_observable("type = backend_info");
	}
}

void MCBackend::print_observables() {
	std::string tmpstr("");
	for(int i = 0; i < _MC_moves; i++) {
		number ratio = (_tries[i] > 0) ? _accepted[i] / (float) _tries[i] : 0;
		//_backend_info += Utils::sformat("  %5.3lf", ratio);
		tmpstr += Utils::sformat("  %5.3lf", ratio);
	}
	_backend_info.insert(0, tmpstr + "  ");

	SimBackend::print_observables();
}
