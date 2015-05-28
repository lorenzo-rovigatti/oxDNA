/*
 * SimManager.cpp
 *
 *  Created on: 03/set/2010
 *      Author: lorenzo
 */

#include <csignal>
#include <unistd.h>
#include "SimManager.h"
#include "../Backends/BackendFactory.h"

void gbl_terminate (int arg) {
	// if the simulation has not started yet, then we make it so pressing ctrl+c twice
	// kills the program no matter what.
	if(!SimManager::started) signal(arg, SIG_DFL);
	OX_LOG(Logger::LOG_INFO, "# Caught SIGNAL %d; setting stop = 1\n", arg);
	SimManager::stop = true;
}

bool SimManager::stop = false;
bool SimManager::started = false;

SimManager::SimManager(int argc, char *argv[]) : _print_energy_every(1000), _restart_step_counter(false) {
	_start_step = _cur_step = _steps = _equilibration_steps = 0;
	_time_scale_manager.state = 0;
	_print_input = 0;
	_pid = getpid();
	_seed = -1;
	_backend = NULL;
	_fix_diffusion_every = 100000;
	_max_steps = -1;
	_time_scale = -1;

	loadInputFile(&_input, argv[1]);
	if(_input.state == ERROR) throw oxDNAException("Caught an error while opening the input file");
	argc -= 2;
	if(argc > 0) addCommandLineArguments(&_input, argc, argv+2);

}

SimManager::~SimManager() {
	// print unread (i.e. unused) keys
	setUnreadKeys(&_input);
	std::string unread;
	for(std::vector<std::string>::iterator it = _input.unread_keys.begin(); it != _input.unread_keys.end(); it++) {
		unread += string("\n\t") + *it;
	}
	if(unread.size() > 0) OX_DEBUG("The following keys found in the input file were not used: %s", unread.c_str());

	cleanInputFile(&_input);
	cleanTimeScale(&_time_scale_manager);

	if(_backend != NULL) {
		int updated = _backend->get_N_updates();
		if(updated > 0) OX_LOG(Logger::LOG_INFO, "Lists updated %d times (every ~%lf steps)", updated, (_cur_step - _start_step) / (double)updated);

		delete _backend;
	}
}

void SimManager::_get_options() {
	getInputBoolAsInt(&_input, "print_input", &_print_input, 0);
	getInputInt(&_input, "print_energy_every", &_print_energy_every, 0);
	getInputInt(&_input, "restart_step_counter", &_restart_step_counter, 0);
	getInputLLInt(&_input, "steps", &_steps, 1);
	getInputLLInt(&_input, "equilibration_steps", &_equilibration_steps, 0);
	if (_equilibration_steps < 0) throw oxDNAException ("Equilibration steps can not be < 0. Aborting");
	if (getInputInt(&_input, "seed", &_seed, 0) == KEY_NOT_FOUND) {
		_seed = time(NULL);
		int rand_seed = 0;
		FILE *f = fopen("/dev/urandom", "rb");
		if (f == NULL) {
			OX_LOG(Logger::LOG_INFO, "Can't open /dev/urandom, using system time as a random seed");
		}
		else {
			if(fread((void *) &rand_seed, sizeof(rand_seed), 1, f) != 0) _seed += rand_seed;
			else OX_LOG(Logger::LOG_INFO, "Can't read from /dev/urandom, using system time as a random seed");
			fclose(f);
		}
	}

	getInputInt(&_input, "fix_diffusion_every", &_fix_diffusion_every, 0);
}

void SimManager::load_options() {
	Logger::instance()->get_settings(_input);
	_get_options();
}

void SimManager::init() {
	OX_LOG(Logger::LOG_INFO, "seeding the RNG with %d", _seed);
	srand48(_seed);

	OX_LOG(Logger::LOG_INFO, "Initializing backend ", _seed);
	_backend = BackendFactory::make_backend(_input);
	_backend->get_settings(_input);

	// here we handle a few SIG* signals;
	signal(SIGTERM, gbl_terminate);
	signal(SIGABRT, gbl_terminate);
	signal(SIGINT, gbl_terminate);
	signal(SIGUSR2, gbl_terminate);

	if(_print_input) {
		char out_name[512];
		sprintf(out_name, "input.%d", _pid);
		printInput(&_input, out_name);
	}
	OX_LOG(Logger::LOG_INFO, "pid: %d", _pid);

	char ts_type[256];
	getInputString(&_input, "time_scale", ts_type, 1);
	_time_scale = TS_LIN;
	if(strcmp(ts_type, "linear") == 0) _time_scale = TS_LIN;
	else if(strcmp(ts_type, "log_lin") == 0) _time_scale = TS_LOG_LIN;
	else throw oxDNAException("Time scale '%s' not supported", ts_type);

	_backend->init();

	if(_restart_step_counter != 0) _start_step = 0;
	else _start_step = _backend->_start_step_from_file;

	// init time_scale_manager
	initTimeScale(&_time_scale_manager, _time_scale);

	int tmp, tmpm;
	getInputInt(&_input, "print_conf_interval", &tmpm, 1);
	setTSInterval(&_time_scale_manager, tmpm);

	if(_time_scale == TS_LOG_LIN) {
		getInputInt(&_input, "print_conf_ppc", &tmp, 1);
		setTSPPC(&_time_scale_manager, tmp);
	}
	// the awkward second argument in the next line makes consistent
	// trajectory files...
	setTSInitialStep(&_time_scale_manager, _start_step + tmpm - (_start_step % tmpm));
	// end

	//_start_step = _start_step + _equilibration_steps;
	_max_steps = _start_step + _steps;
}

void SimManager::run() {
	SimManager::started = true;
	// equilibration loop
	if(_equilibration_steps > 0) {
		OX_LOG(Logger::LOG_INFO, "Equilibrating...");
		for(llint step = 0; step < _equilibration_steps && !SimManager::stop; step++) {
			_backend->sim_step(step);
			if (step > 1 && step % _fix_diffusion_every == 0) _backend->fix_diffusion();
		}
		OX_LOG(Logger::LOG_INFO, "Equilibration done");
	}

	// main loop
	for(_cur_step = _start_step; _cur_step < _max_steps && !SimManager::stop; _cur_step++) {
		if(_cur_step == _time_scale_manager.next_step) {
			if(_cur_step > _start_step) _backend->print_conf(_cur_step);
			setTSNextStep(&_time_scale_manager);
		}
		if (_cur_step > 1 && _cur_step % _fix_diffusion_every == 0) _backend->fix_diffusion();

		_backend->print_observables(_cur_step);
		_backend->sim_step(_cur_step);
	}
	// this is in case _cur_step, after being increased by 1 before exiting the loop,
	// has become a multiple of print_conf_every
	if (_cur_step > 1 && _cur_step % _fix_diffusion_every == 0) _backend->fix_diffusion();

	if(_cur_step == _time_scale_manager.next_step) {
		if(_cur_step > _start_step) _backend->print_conf(_cur_step);
		setTSNextStep(&_time_scale_manager);
	}
	_backend->print_observables(_cur_step);

	// prints the last configuration
	_backend->print_conf(_cur_step, false, true);

	TimingManager::instance()->print(_cur_step - _start_step);
}

