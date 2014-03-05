/*
 * CSDManager.cpp
 *
 *  Created on: 26/set/2011
 *      Author: lorenzo
 */

#include <string>

#include "CSDManager.h"

CSDManager::CSDManager(IOManager &IO, const char *input_file) : SimManager(IO, input_file), _N_confs(0), _confs_to_skip(0) {
	_csdbackend = new CSDBackend(&IO);
	// we do this in order to not to break backend's initialization
	_backend = _csdbackend;
}

CSDManager::~CSDManager() {
	delete[] _csd;
}

void CSDManager::_get_options() {
	SimManager::_get_options();

	getInputString(&_input, "csd_output", _csd_output, 1);
	getInputInt(&_input, "confs_to_skip", &_confs_to_skip, 0);
}

void CSDManager::_load_backend(int confs_to_skip) {
	delete _csdbackend;
	_csdbackend = new CSDBackend(&_IO);
	_csdbackend->get_settings(_input);

	ifstream conf_input(_conf_file);
	_csdbackend->set_confs_to_skip(confs_to_skip);
	_csdbackend->init(conf_input);
	conf_input.close();

	_backend = _csdbackend;
}

void CSDManager::init() {
	SimManager::init();

	_N_cyls = _csdbackend->_N_cyls;
	_csd = new int[_N_cyls+1];
	for(int i = 0; i <= _N_cyls; i++) _csd[i] = 0;

	// calculate the number of configurations in the trajectory file
	string line;
	ifstream conf_input(_conf_file);
	while(getline(conf_input, line)) if(line[0] == 'E') _N_confs++;
	conf_input.close();

	_N_confs -= _confs_to_skip;
	_IO.log(_IO.LOG_INFO, "Number of configurations: %d", _N_confs);
}

void CSDManager::run() {
	for(int i = 0; i < _N_confs && !SimManager::stop; i++) {
		_load_backend(_confs_to_skip + i);
		_csdbackend->set_csd();

		for(int j = 0; j <= _N_cyls; j++) {
			_csd[j] += _csdbackend->_csd[j];
		}
	}

	ofstream out(_csd_output);
	for(int i = 1; i <= _N_cyls; i++) {
		if(_csd[i] > 0) out << i << " " << (double)_csd[i]/(double)_N_confs << endl;
	}
	out.close();
}
