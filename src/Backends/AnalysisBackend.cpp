/*
 * AnalysisBackend.cpp
 *
 *  Created on: Feb 13, 2013
 *      Author: rovigatti
 */

#include <sstream>

#include "AnalysisBackend.h"
#include "../Interactions/InteractionFactory.h"
#include "../Observables/ObservableOutput.h"
#include "../Lists/ListFactory.h"
#include "../Boxes/BoxFactory.h"
#include "../PluginManagement/PluginManager.h"

AnalysisBackend::AnalysisBackend() :
				SimBackend(),
				_done(false),
				_n_conf(0) {
	_enable_fix_diffusion = 0;
}

AnalysisBackend::~AnalysisBackend() {

}

void AnalysisBackend::get_settings(input_file &inp) {
	_config_info->sim_input = &inp;

	// initialise the plugin manager with the input file
	PluginManager::instance()->init(inp);

	_interaction = InteractionFactory::make_interaction(inp);
	_interaction->get_settings(inp);

	_box = BoxFactory::make_box(inp);
	_box->get_settings(inp);

	_lists = ListFactory::make_list(inp, _particles, _box.get());
	_lists->get_settings(inp);

	// initialise the timer
	_mytimer = TimingManager::instance()->new_timer(std::string("AnalysisBackend"));

	getInputInt(&inp, "analysis_confs_to_skip", &_confs_to_skip, 0);

	getInputString(&inp, "trajectory_file", _conf_filename, 1);

	getInputDouble(&inp, "max_io", &_max_io, 0);

	getInputBool(&inp, "binary_initial_conf", &_initial_conf_is_binary, 0);
	if(_initial_conf_is_binary) {
		OX_LOG(Logger::LOG_INFO, "Reading binary configuration");
	}

	int tmp;
	int val = getInputBoolAsInt(&inp, "external_forces", &tmp, 0);
	if(val == KEY_FOUND) {
		_external_forces = (tmp != 0);
		if(_external_forces) {
			getInputString(&inp, "external_forces_file", _external_filename, 1);
		}
	}
	else if(val == KEY_INVALID) throw oxDNAException("external_forces must be either 0 (false, no) or 1 (true, yes)");

	char raw_T[256];
	getInputString(&inp, "T", raw_T, 1);
	_T = Utils::get_temperature(raw_T);

	// here we fill the _obs_outputs vector
	int i = 1;
	bool found = true;
	while(found) {
		stringstream ss;
		ss << "analysis_data_output_" << i;
		string obs_string;
		if(getInputString(&inp, ss.str().c_str(), obs_string, 0) == KEY_FOUND) {
			ObservableOutputPtr new_obs_out = std::make_shared<ObservableOutput>(obs_string);
			add_output(new_obs_out);
		}
		else found = false;

		i++;
	}
}

void AnalysisBackend::init() {
	SimBackend::init();
}

const FlattenedConfigInfo &AnalysisBackend::flattened_conf() {
	_flattened_conf.update(_read_conf_step, particles());
	return _flattened_conf;
}

bool AnalysisBackend::read_next_configuration(bool binary) {
	_done = !SimBackend::read_next_configuration(binary);
	return !_done;
}

void AnalysisBackend::analyse() {
	_mytimer->resume();

	if(_n_conf % 100 == 0 && _n_conf > 0) {
		OX_LOG(Logger::LOG_INFO, "Analysed %d configurations", _n_conf);
	}
	SimBackend::print_observables(_read_conf_step);

	if(read_next_configuration(_initial_conf_is_binary)) {
		_n_conf++;

		for(auto p : _particles) {
			_lists->single_update(p);
		}
		_lists->global_update();
	}

	_mytimer->pause();
}
