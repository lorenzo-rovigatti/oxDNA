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

AnalysisBackend::AnalysisBackend() : SimBackend<double>(), _done(false), _n_conf(0) {
	_enable_fix_diffusion = 0;
}

AnalysisBackend::~AnalysisBackend() {

}

void AnalysisBackend::get_settings(input_file &inp) {
	// initialise the plugin manager with the input file
	PluginManager *pm = PluginManager::instance();
	pm->init(inp);

	_interaction = InteractionFactory::make_interaction<double>(inp);
	_interaction->get_settings(inp);
	
	_box = BoxFactory::make_box<double>(inp);
	_box->get_settings(inp);

	_lists = ListFactory::make_list<double>(inp, _N, _box);
	_lists->get_settings(inp);

	// initialise the timer
	_mytimer = TimingManager::instance()->new_timer(std::string("AnalysisBackend"));

	getInputInt(&inp, "analysis_confs_to_skip", &_confs_to_skip, 0);

	getInputString(&inp, "trajectory_file", this->_conf_filename, 1); 

	getInputDouble(&inp, "max_io", &this->_max_io, 0);

	getInputBool(&inp, "binary_initial_conf", &_initial_conf_is_binary, 0);
	if (_initial_conf_is_binary){
  	OX_LOG(Logger::LOG_INFO, "Reading binary configuration");
  }

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

void AnalysisBackend::init() {
	SimBackend<double>::init();
}

void AnalysisBackend::analyse() {
	_mytimer->resume();

	if(_n_conf % 100 == 0 && _n_conf > 0) OX_LOG(Logger::LOG_INFO, "Analysed %d configurations", _n_conf);
	SimBackend<double>::print_observables(_read_conf_step);

	if(!_read_next_configuration(_initial_conf_is_binary)) _done = true;
	else _n_conf++;

	for(int i = 0; i < this->_N; i++) this->_lists->single_update(this->_particles[i]);
	this->_lists->global_update();

	_mytimer->pause();
}
