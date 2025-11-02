/*
 * BaseObservable.h
 *
 *  Created on: Mar 21, 2020
 *      Author: lorenzo
 */

#ifndef OXPY_BUILD_DOCS_BINDINGS_INCLUDES_BASEOBSERVABLE_H_
#define OXPY_BUILD_DOCS_BINDINGS_INCLUDES_BASEOBSERVABLE_H_

#include "../../python_defs.h"

#include "HBList.h"

#include <Observables/BaseObservable.h>

// trampoline class for BaseObservable
class PyBaseObservable : public BaseObservable {
public:
	using BaseObservable::BaseObservable;

	std::string get_output_string(llint curr_step) override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				std::string,
				BaseObservable,
				get_output_string,
				curr_step
		);

		return std::string();
	}
};

void export_BaseObservable(py::module &m) {
	py::module sub_m = m.def_submodule("observables");

	py::class_<BaseObservable, PyBaseObservable, std::shared_ptr<BaseObservable>> obs(sub_m, "BaseObservable", R"pbdoc(
		The interface class for observables.
	)pbdoc");

	obs.def(py::init<>(), R"pbdoc(
        The default constructor takes no parameters.
	)pbdoc");

	obs.def("get_settings", &BaseObservable::get_settings, py::arg("my_inp"), py::arg("sim_inp"), R"pbdoc(
		Get the settings required to initialise the observable from the simulation and its own input files.

		Parameters
		---------- 
		my_inp : :class:`input_file`
			The input file of the observable.
		sim_inp : :class:`input_file`
			The general input file of the simulation.
	)pbdoc");

	obs.def("init", &BaseObservable::init, R"pbdoc(
		Initialise the observable.
	)pbdoc");

	obs.def_property_readonly("config_info", [](BaseObservable const& obs) { return CONFIG_INFO->instance(); }, R"pbdoc(
        The simulation's :class:`ConfigInfo` object, which can be used from within :meth:`get_output_string` to 
        access the current configuration's details.
	)pbdoc");

	obs.def_property("id", &BaseObservable::get_id, &BaseObservable::set_id, "The id of the observable.");

	obs.def("get_output_string", &BaseObservable::get_output_string, py::arg("curr_step"), R"pbdoc(
		Compute the quantity/quantities of interest and returns the output string.

        Parameters
        ---------- 
        curr_step : int
            The current simulation step.

        Returns
        -------
        str
            The output of the observable.
	)pbdoc");

	export_HBList(sub_m);
}

#endif /* OXPY_BUILD_DOCS_BINDINGS_INCLUDES_BASEOBSERVABLE_H_ */
