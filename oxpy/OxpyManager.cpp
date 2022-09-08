/*2
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
input_file _input_file_from_filename(std::string input_filename) {
	input_file input;
	input.init_from_filename(input_filename.c_str());

	return input;
}

OxpyManager::OxpyManager(std::string input_filename) :
				OxpyManager(_input_file_from_filename(input_filename)) {

}

OxpyManager::OxpyManager(input_file input) :
				SimManager(input) {
SimManager::load_options();
SimManager::init();
}

OxpyManager::~OxpyManager() {

}

std::shared_ptr<ConfigInfo> OxpyManager::config_info() {
	return CONFIG_INFO;
}

number OxpyManager::system_energy() {
	return CONFIG_INFO->interaction->get_system_energy(CONFIG_INFO->particles(), CONFIG_INFO->lists);
}

void OxpyManager::update_temperature(number new_T) {
	CONFIG_INFO->update_temperature(new_T);
}

void OxpyManager::print_configuration(bool also_last) {
	// prints the trajectory configuration
	_backend->print_conf();
	// prints the last configuration
	_backend->print_conf(false, true);
}

void OxpyManager::add_output(std::string filename, llint print_every, std::vector<ObservablePtr> observables) {
	std::string output_string = Utils::sformat("{\n\tname = %s\n\tprint_every = %lld\n}\n", filename.c_str(), print_every);

	ObservableOutputPtr output = std::make_shared<ObservableOutput>(output_string);
	for(auto obs : observables) {
		output->add_observable(obs);
	}
	output->init();

	_backend->add_output(output);
}

void OxpyManager::remove_output(std::string filename) {
	_backend->remove_output(filename);
}

void OxpyManager::update_CPU_data_structures() {
	_backend->apply_simulation_data_changes();
}

void OxpyManager::run(llint steps, bool print_output) {
	_backend->apply_changes_to_simulation_data();

	for(llint i = 0; i < steps && !SimManager::stop; i++, _steps_run++) {
		if(_backend->current_step() == _time_scale_manager.next_step) {
			if(print_output && i > 0) {
				_backend->print_conf();
			}
			setTSNextStep(&_time_scale_manager);
		}

		if(i > 0 && i % _fix_diffusion_every == 0) {
			_backend->fix_diffusion();
		}

		_backend->update_observables_data();
		if(print_output) {
			_backend->print_observables();
		}

		_backend->sim_step();
		_backend->increment_current_step();
	}

	_backend->apply_simulation_data_changes();
}

void OxpyManager::print_timings() {
	TimingManager::instance()->print(_steps_run);
}

void export_SimManager(py::module &m) {
	pybind11::class_<SimManager, std::shared_ptr<SimManager>> manager(m, "SimManager");

	manager.def("run_complete", &SimManager::run, R"pbdoc(
        Run the simulation till completion (that is, till the simulation has run number of steps specified in the input file).
	)pbdoc");
}

void export_OxpyManager(py::module &m) {
	py::class_<OxpyManager, SimManager, std::shared_ptr<OxpyManager>> manager(m, "OxpyManager", R"pbdoc(
		The object in charge of initialising and running a simulation. 
	)pbdoc");

	manager.def(py::init<std::string>(), py::arg("input_filename"), R"pbdoc(
The default constructor takes the input filename as its only parameter.

Parameters
----------
input_filename : str
    The name of the input file that will be used to run the simulation
	)pbdoc");

	manager.def(py::init<input_file>(), py::arg("input"), R"pbdoc(
This constructor takes an :class:`InputFile` as its only argument.

Parameters
----------
input: :class:`InputFile`
    The object storing the simulations' options.
	)pbdoc");

	manager.def("load_options", [](OxpyManager &manager) { OX_LOG(Logger::LOG_WARNING, "OxpyManager.load_options() is deprecated and will be removed in future releases since it is automatically called by the class' constructor."); }, R"pbdoc(
		.. deprecated:: 3.2.2

		This method does nothing and it exists only for backward compatibility reasons.
	)pbdoc");

	manager.def("init", [](OxpyManager &manager) { OX_LOG(Logger::LOG_WARNING, "OxpyManager.init() is deprecated and will be removed in future releases since it is automatically called by the class' constructor."); }, R"pbdoc(
		.. deprecated:: 3.2.2

		This method does nothing and it exists only for backward compatibility reasons.
	)pbdoc");

	manager.def("print_configuration", &OxpyManager::print_configuration, py::arg("also_last") = true, R"pbdoc(
        Append the current configuration to the trajectory file and, by default, print it in the "lastconf_file".

        Parameters
        ----------
            also_last : bool
                If True (default value) prints the current configuration also in the "lastconf_file" file.  
	)pbdoc");

	manager.def("config_info", &OxpyManager::config_info, R"pbdoc(
		Return the simulation's :class:`ConfigInfo`.

        Returns
        -------
        :class:`ConfigInfo`
            The object that stores all the simulation's details (particles, interaction, `etc`).
	)pbdoc");

	manager.def("system_energy", &OxpyManager::system_energy, R"pbdoc(
		Compute the potential energy of the system in its current state.

        Returns
        -------
        float
            The system's potential energy.
	)pbdoc");

	manager.def("update_temperature", &OxpyManager::update_temperature, pybind11::arg("new_T"), R"pbdoc(
        Change the temperature at which the simulation is run.

        Parameters
        ----------
            new_T : float
                The new temperature.
	)pbdoc");

	manager.def("add_output", &OxpyManager::add_output, pybind11::arg("filename"), pybind11::arg("print_every"), pybind11::arg("observables"), R"pbdoc(
        Add an output file to the simulation, to be populated with the given observables. Added outputs can be subsequently removed with :meth:`remove_output`.

        Parameters
        ----------
            filename : str
				The name of the file where the output will be printed.
            print_every : int
                The frequency (in simulation time steps) with which the output should be computed and printed. 
            observables : List(:class:`BaseObservable`)
                A list of observables that will be used to populate the output.
	)pbdoc");

	manager.def("remove_output", &OxpyManager::remove_output, pybind11::arg("filename"), R"pbdoc(
		Remove an output file from the simulation. Once an output is removed, the simulation will stop updating it.

		Parameters
		----------
			filename : str
				The name of the output file, which should be the same name used to add it (for instance *via* :meth:`add_output`).
	)pbdoc");

	manager.def("update_CPU_data_structures", &OxpyManager::update_CPU_data_structures, R"pbdoc(
		Update the CPU data structures. Useful only when running CUDA simulations.

		This method copies simulation data from the GPU to the CPU, and as such calling it too often may severely decrease performance.
		It is automatically invoked by meth:`run`, and therefore it makes sense to call it only in specific cases (*e.g.* to access simulation data from callbacks).
	)pbdoc");

	manager.def("run", &OxpyManager::run, pybind11::arg("steps"), pybind11::arg("print_output") = true, R"pbdoc(
		Run the simulation for the given number of steps. The second argument controls whether the simulations output (configurations and observables) should be printed or not.

		Parameters
		----------
			steps : int
				The number of steps to be run.
			print_output : bool
				If True (the default value) the simulation output will be printed.
	)pbdoc");

	manager.def("steps_run", &OxpyManager::steps_run, R"pbdoc(
		Return the number of steps run using this manager.
		
		Returns
		-------
			int
				The number of steps run.
	)pbdoc");

	manager.def("print_timings", &OxpyManager::print_timings, R"pbdoc(
		Print the timings taking into account the number of steps run by the manager.
	)pbdoc");

	manager.def_property_readonly("current_step", [](OxpyManager &manager) { return manager.config_info()->curr_step; }, "The current time step.");
}
