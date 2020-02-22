/*
 * OxpyManager.cpp
 *
 *  Created on: Sep 13, 2019
 *      Author: lorenzo
 */

#include "OxpyManager.h"

/**
 * Converts the argument to a vector of char * to be consumed by the SimManager constructor. The first element
 * of the new vector is initialised to "oxpy" because SimManager expectes an argv-like vector, with the first element
 * containing the name of the program.
 *
 * We need this because pybind11 does not support binding methods that take double pointers.
 *
 * @param to_convert
 * @return
 */
char **_convert_to_SimManager_argv(std::vector<std::string> to_convert) {
	static char exe_name[10];
	sprintf(exe_name, "oxpy");
	std::vector<char *> argv;
	argv.reserve(to_convert.size() + 1);
	argv.push_back(exe_name);
	for(auto &str : to_convert) {
		argv.push_back(const_cast<char *>(str.c_str()));
	}

	return argv.data();
}

OxpyManager::OxpyManager(std::vector<std::string> args) :
				SimManager(args.size() + 1, _convert_to_SimManager_argv(args)) {

}

OxpyManager::~OxpyManager() {

}

std::shared_ptr<ConfigInfo> OxpyManager::config_info() {
	return CONFIG_INFO;
}

void OxpyManager::run(llint steps, bool print_output) {
	if(_cur_step < _start_step) {
		_cur_step = _start_step;
	}

	_backend->apply_changes_to_simulation_data();

	for(llint i = 0; i < steps && !SimManager::stop; i++, _cur_step++) {
		if(_cur_step == _time_scale_manager.next_step) {
			if(print_output && _cur_step > _start_step) {
				_backend->print_conf(_cur_step);
			}
			setTSNextStep(&_time_scale_manager);
		}

		if(_cur_step > 0 && _cur_step % _fix_diffusion_every == 0) {
			_backend->fix_diffusion();
		}

		if(print_output) {
			_backend->print_observables(_cur_step);
		}

		_backend->sim_step(_cur_step);
	}

	_backend->apply_simulation_data_changes();
}

void export_SimManager(py::module &m) {
	pybind11::class_<SimManager, std::shared_ptr<SimManager>> manager(m, "SimManager");

	manager
		.def("load_options", &SimManager::load_options)
		.def("init", &SimManager::init)
		.def("run_complete", &SimManager::run);
}

void export_OxpyManager(py::module &m) {
	pybind11::class_<OxpyManager, SimManager, std::shared_ptr<OxpyManager>> manager(m, "OxpyManager");

	manager
		.def(py::init<std::vector<std::string>>())
		.def("config_info", &OxpyManager::config_info)
		.def("run", &OxpyManager::run, pybind11::arg("steps"), pybind11::arg("print_output") = true);
}
