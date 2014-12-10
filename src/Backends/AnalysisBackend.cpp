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
#include "../PluginManagement/PluginManager.h"

AnalysisBackend::AnalysisBackend() : SimBackend<double>(), _done(false), _n_conf(0) {
	_enable_fix_diffusion = 0;
}

AnalysisBackend::~AnalysisBackend() {
	PluginManager::clear();
}

void AnalysisBackend::get_settings(input_file &inp) {
	// initialise the plugin manager with the input file
	PluginManager *pm = PluginManager::instance();
	pm->init(inp);

	_interaction = InteractionFactory::make_interaction<double>(inp);
	_interaction->get_settings(inp);

	_lists = ListFactory::make_list<double>(inp, _N, _box_side);
        _lists->get_settings(inp);

	getInputInt(&inp, "analysis_confs_to_skip", &_confs_to_skip, 0);

	int tmp;
	int val = getInputBoolAsInt(&inp, "external_forces", &tmp, 0);
	if(val == KEY_FOUND) {
		_external_forces = (tmp != 0);
		if (_external_forces) {
			getInputString (&inp, "external_forces_file", _external_filename, 1);
		}
	}
	else if(val == KEY_INVALID) throw oxDNAException("external_forces must be either 0 (false, no) or 1 (true, yes)");

	char raw_T[256];
	getInputString(&inp, "T", raw_T, 1);
	char deg;
	double tmp_T;
	int res = sscanf(raw_T, "%lf %c", &tmp_T, &deg);
	if(res == 2) {
		deg = tolower(deg);
		switch(deg) {
		case 'c':
			_T = (tmp_T + 273.15) * 0.1 / 300.; // convert to kelvin and then to simulation units
			OX_LOG(Logger::LOG_INFO, "Converting temperature from Celsius (%lf C°) to simulation units (%lf)", tmp_T, _T);
			break;
		case 'k':
			_T = tmp_T * 0.1 / 300.; // convert to simulation units
			OX_LOG(Logger::LOG_INFO, "Converting temperature from Kelvin (%lf K) to simulation units (%lf)", tmp_T, _T);
			break;
		default:
			throw oxDNAException("Unrecognizable temperature '%s'", raw_T);
			/* no break */
		}
	}
	else _T = tmp_T;

	// here we fill the _obs_outputs vector
	int i = 1;
	bool found = true;
	while(found) {
		stringstream ss;
		ss << "analysis_data_output_" << i;
		string obs_string;
		if(getInputString(&inp, ss.str().c_str(), obs_string, 0) == KEY_FOUND) {
			ObservableOutput<double> *new_obs_out = new ObservableOutput<double>(obs_string, inp);
			_obs_outputs.push_back(new_obs_out);
		}
		else found = false;

		i++;
	}
}

void AnalysisBackend::init(char traj_path[256]) {
	SimBackend<double>::init(traj_path);
	// this is to avoid the "print timings" at the end
	destroy_timer(&_timer);
}

void AnalysisBackend::analyse() {
	if(_n_conf % 100 == 0 && _n_conf > 0) OX_LOG(Logger::LOG_INFO, "Analysed %d configurations", _n_conf);
	SimBackend<double>::print_observables(_read_conf_step);
	if(!_read_next_configuration()) _done = true;
	else _n_conf++;
}
