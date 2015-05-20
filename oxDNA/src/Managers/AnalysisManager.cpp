/*
 * AnalysisManager.cpp
 *
 *  Created on: Feb 13, 2013
 *      Author: rovigatti
 */

#include "AnalysisManager.h"

AnalysisManager::AnalysisManager(int argc, char *argv[]) {
	loadInputFile(&_input, argv[1]);
	if(_input.state == ERROR) throw oxDNAException("Caught an error while opening the input file");
	argc -= 2;
	if(argc > 0) addCommandLineArguments(&_input, argc, argv+2);

	_backend = new AnalysisBackend();
}

AnalysisManager::~AnalysisManager() {
	cleanInputFile(&_input);
	delete _backend;
}

void AnalysisManager::load_options() {
	Logger::instance()->get_settings(_input);

	_backend->get_settings(_input);
}

void AnalysisManager::init() {
	_backend->init();
}

void AnalysisManager::analysis() {
	while(!_backend->done()) _backend->analyse();
}
