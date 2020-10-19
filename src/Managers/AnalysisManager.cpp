/*
 * AnalysisManager.cpp
 *
 *  Created on: Feb 13, 2013
 *      Author: rovigatti
 */

#include "AnalysisManager.h"

#include "../Utilities/oxDNAException.h"

AnalysisManager::AnalysisManager(input_file input) :
				_input(input) {
	_backend = std::make_shared<AnalysisBackend>();
}

AnalysisManager::~AnalysisManager() {

}

void AnalysisManager::load_options() {
	Logger::instance()->get_settings(_input);

	// seed;
	int seed;
	if(getInputInt(&_input, "seed", &seed, 0) == KEY_NOT_FOUND) {
		seed = time(NULL);
		int rand_seed = 0;
		FILE *f = fopen("/dev/urandom", "rb");
		if(f == NULL) {
			OX_LOG(Logger::LOG_INFO, "Can't open /dev/urandom, using system time as a random seed");
		}
		else {
			if(fread((void *) &rand_seed, sizeof(rand_seed), 1, f) != 0) seed += rand_seed;
			else OX_LOG(Logger::LOG_INFO, "Can't read from /dev/urandom, using system time as a random seed");
			fclose(f);
		}
	}
	OX_LOG(Logger::LOG_INFO, "Setting the random number generator with seed = %d", seed);
	srand48((long int) seed);

	_backend->get_settings(_input);
}

void AnalysisManager::init() {
	_backend->init();
}

void AnalysisManager::analysis() {
	while(!_backend->done()) {
		_backend->analyse();
	}
}
