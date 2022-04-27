/*
 * AnalysisBackend.h
 *
 *  Created on: 17 feb 2021
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_ANALYSISBACKEND_H_
#define OXPY_BINDINGS_INCLUDES_ANALYSISBACKEND_H_

#include "../python_defs.h"

#include <Backends/AnalysisBackend.h>

std::shared_ptr<AnalysisBackend> backend_builder(input_file &input) {
	std::shared_ptr<AnalysisBackend> backend = std::make_shared<AnalysisBackend>();
	backend->get_settings(input);
	backend->init();

	return backend;
}

void export_AnalysisBackend(py::module &m) {
	py::class_<AnalysisBackend, std::shared_ptr<AnalysisBackend>> backend(m, "AnalysisBackend", R"pbdoc(
		An object that can be used to analyse oxDNA trajectory files.
	)pbdoc");

	// here we provide a "builder" that creates the object and initialises it, so that the user doesn't have to explicitly do it
	backend.def(py::init(&backend_builder), R"pbdoc(
        The constructor takes an :py:class:`~oxpy.core.InputFile` that specifies the options needed to build the backend.
	)pbdoc");

	backend.def("read_next_configuration", &AnalysisBackend::read_next_configuration, py::arg("binary")=false, R"pbdoc(
		Load up the next configuration from the trajectory file.
	)pbdoc");

	backend.def("analyse", &AnalysisBackend::analyse, R"pbdoc(
		Analyse the current configuration using the observables set through the input file.
	)pbdoc");

	backend.def("done", &AnalysisBackend::done, R"pbdoc(
        True if the backend has reached the last configuration in the trajectory file, False otherwise.

        Returns
        -------
        bool
            The completion status of the analysis.
	)pbdoc");

	backend.def_property_readonly("particles", &AnalysisBackend::particles, R"pbdoc(
		The list of all the particles.
	)pbdoc");

	backend.def_property_readonly("flattened_conf", &AnalysisBackend::flattened_conf, R"pbdoc(
		 A flattened (that is, composed by arrays rather than class instances) view of the current configuration, stored in a :py:class:`~oxpy.core.FlattenedConfigInfo` object.
	)pbdoc");

	backend.def("config_info", &AnalysisBackend::config_info, R"pbdoc(
		Return the simulation's :class:`ConfigInfo`.

		Returns
		-------
		:class:`ConfigInfo`
			The object that stores all the simulation's details (particles, interaction, `etc`).
	)pbdoc");

	backend.def_property_readonly("conf_step", &AnalysisBackend::get_conf_step, "The time step at which the current configuration was printed.");
}

#endif /* OXPY_BINDINGS_INCLUDES_ANALYSISBACKEND_H_ */
