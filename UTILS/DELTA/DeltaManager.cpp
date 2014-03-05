/*
 * DeltaManager.cpp
 *
 *  Created on: 20/set/2011
 *      Author: lorenzo
 */

#include "DeltaManager.h"

DeltaManager::DeltaManager(IOManager &IO, const char *input_file) : SimManager(IO, input_file), _append(0) {
	_dbackend = new DeltaBackend(&IO);
	// we do this in order to not to break backend's initialization
	_backend = _dbackend;
}

DeltaManager::~DeltaManager() {

}

void DeltaManager::_get_options() {
	SimManager::_get_options();

	getInputInt(&_input, "change_every", &_change_every, 1);
	getInputInt(&_input, "every_change_skip", &_every_change_skip, 1);
	getInputInt(&_input, "confs_to_skip", &_confs_to_skip, 1);
	getInputInt(&_input, "delta_append", &_append, 0);
	if(getInputDouble(&_input, "delta_threshold", &_threshold, 0) == KEY_NOT_FOUND) _threshold = -1e-8;
	getInputString(&_input, "E_output", _E_output, 1);
}

void DeltaManager::run() {
	int found;
	FILE *E_out;
	if(_append) E_out = fopen(_E_output, "a");
	else E_out = fopen(_E_output, "w");

	fprintf(E_out, "%lf\n", _dbackend->get_trial_volume());
	_cur_step = 1;
	int n_conf = 0;

	while(!SimManager::stop) {
		found = 0;
		while(found < _change_every && !SimManager::stop) {
			double E = _dbackend->insertion_energy();
			if (E < _threshold) {
				found++;
				fprintf(E_out, "%lld %e\n", _cur_step, E);
				//fprintf(stdout, "%lld %e\n", _cur_step, E);
				fflush(E_out);
			}

			_cur_step++;
		}

		if(!SimManager::stop) {
			n_conf++;

			delete _dbackend;
			_dbackend = new DeltaBackend(&_IO);
			_dbackend->get_settings(_input);
			_dbackend->set_confs_to_skip(_confs_to_skip + n_conf*_every_change_skip);
			ifstream conf_input(_conf_file);
			_dbackend->init(conf_input);
			conf_input.close();
			
			_backend = _dbackend;
			fprintf(E_out, "%lf\n", _dbackend->get_trial_volume());
		}
	}

	_dbackend->print_conf(_cur_step, false, true);

	fclose(E_out);
}
