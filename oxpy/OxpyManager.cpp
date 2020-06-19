/*
 * OxpyManager.cpp
 *
 *  Created on: Sep 13, 2019
 *      Author: lorenzo
 */

#include "OxpyManager.h"

/**
 * First converts the argument to an array of char * and then builds the input_file with it. The first element of the new array is
 * initialised to "oxpy" because SimManager expectes an argv-like vector, with the first element containing the name of the program.
 *
 * @param to_convert
 * @return
 */
input_file *input_file_from_args(std::vector<std::string> to_convert) {
	static char exe_name[10];
	sprintf(exe_name, "oxpy");
	std::vector<char *> argv;
	argv.reserve(to_convert.size() + 1);
	argv.push_back(exe_name);
	for(auto &str : to_convert) {
		argv.push_back(const_cast<char *>(str.c_str()));
	}

	input_file *input = new input_file();
	input->init_from_command_line_args(argv.size(), argv.data());

	return input;
}

OxpyManager::OxpyManager(std::vector<std::string> args) :
				SimManager(*input_file_from_args(args)) {

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

	manager.def("load_options", &SimManager::load_options, R"pbdoc(
		Load the options from the input file.
	)pbdoc");

	manager.def("init", &SimManager::init, R"pbdoc(
		Initialise the simulation manager.
	)pbdoc");

	manager.def("run_complete", &SimManager::run, R"pbdoc(
        Run the simulation till completion (that is, till the simulation has run number of steps specified in the input file).
	)pbdoc");
}

void export_OxpyManager(py::module &m) {
	py::class_<OxpyManager, SimManager, std::shared_ptr<OxpyManager>> manager(m, "OxpyManager", R"pbdoc(
		The object in charge of initialising and running a simulation. 
	)pbdoc");

	manager.def(py::init<std::vector<std::string>>(), py::arg("args"), R"pbdoc(
		The default constructor takes a list of strings as its only parameter.

        Parameters
        ----------
        args: List(str)
            The list or arguments to be passed to the manager. The first element should be the name of the input file,
            while the others should be key=value pairs.
	)pbdoc");

	manager.def("config_info", &OxpyManager::config_info, R"pbdoc(
		Return the simulation's :class:`ConfigInfo`.

        Returns
        -------
        :class:`ConfigInfo`
            The object that stores all the simulation's details (particles, interaction, `etc`).
	)pbdoc");

	manager.def("run", &OxpyManager::run, pybind11::arg("steps"), pybind11::arg("print_output") = true, R"pbdoc(
		Run the simulation for the given number of steps. The second argument controls whether the simulations output (configurations and observables) should be printed or not.

        Parameters
        ----------
            steps: int
                The number of steps to be run.
            print_output: bool
                If True (the default value) the simulation output will be printed.
	)pbdoc");
}
