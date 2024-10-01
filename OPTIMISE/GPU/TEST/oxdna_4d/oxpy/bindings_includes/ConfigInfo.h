/*
 * ConfigInfo.h
 *
 *  Created on: Mar 21, 2020
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_CONFIGINFO_H_
#define OXPY_BINDINGS_INCLUDES_CONFIGINFO_H_

#include <Utilities/ConfigInfo.h>
#include "../python_defs.h"

void export_ConfigInfo(py::module &m) {
	py::class_<ConfigInfo, std::shared_ptr<ConfigInfo>> conf_info(m, "ConfigInfo", R"pbdoc(
		This singleton object stores all the details of the simulation (particles, neighbour lists, input file, interaction, external forces) 
	)pbdoc");

	conf_info.def("notify", &ConfigInfo::notify, py::arg("event"), R"pbdoc(
		Notify the triggering of an event. Any callback associated to the event will be invoked.

		Parameters
		----------
		event: str
			The triggered event.
	)pbdoc");

	conf_info.def("subscribe", &ConfigInfo::subscribe, py::arg("event"), py::arg("callback"), R"pbdoc(
		Assign a callback to the given event. 

		The callback will be invoked every time the event is triggered.

		Parameters
		----------
		event: str
			The event associated to the callback.
		callback: callable
			A callable that takes no parameters.

	)pbdoc");

	conf_info.def("N", &ConfigInfo::N, R"pbdoc(
        Return the current number of particles.

        Returns
        -------
        int
            The number of particles in the simulation box.
	)pbdoc");

	conf_info.def("particles", &ConfigInfo::particles, py::return_value_policy::reference, R"pbdoc(
		Return a list of all the particles.

        Returns
        -------
        List(:py:class:`BaseParticle`)
            A list containing all the particles in the simulation box.
	)pbdoc");

	conf_info.def("molecules", &ConfigInfo::molecules, py::return_value_policy::reference, R"pbdoc(
		Return a list of all the molecules.

		Returns
		-------
		List(:py:class:`Molecule`)
			A list containing all the molecules.
	)pbdoc");

	conf_info.def_readwrite("interaction", &ConfigInfo::interaction, R"pbdoc(
		The simulation's :py:class:`IBaseInteraction` object.
	)pbdoc");

	conf_info.def_property_readonly("flattened_conf", &ConfigInfo::flattened_conf, R"pbdoc(
		A flattened (that is, composed by arrays rather than class instances) view of the current configuration, stored in a :py:class:`FlattenedConfigInfo` object.
	)pbdoc");

	conf_info.def_readonly("forces", &ConfigInfo::forces, "The external forces.");

	conf_info.def("get_force_by_id", &ConfigInfo::get_force_by_id, R"pbdoc(
		Return the force with the given id, which can be set in the external forces file with `id = FORCE_ID`. 

		Returns
		-------
		:py:class:`~oxpy.core.forces.BaseForce`
			The external force with the given id, or `None` if the id does not correspond to any external force.
	)pbdoc");

	conf_info.def_readonly("observables", &ConfigInfo::observables, "The simulation observables.");

	conf_info.def("get_observable_by_id", &ConfigInfo::get_observable_by_id, R"pbdoc(
		Return the observable with the given id, which can be set in the portion of the input file relative to the observable `id = OBSERVABLE_ID`. 

		Returns
		-------
		:py:class:`~oxpy.core.observables.BaseObservable`
			The observable with the given id, or `None` if the id does not correspond to any observable.
	)pbdoc");

	conf_info.def_property_readonly("box_sides", [](ConfigInfo &c) { return c.box->box_sides(); }, R"pbdoc(
		The length of the edges of the simulation box.
	)pbdoc");

	conf_info.def_property_readonly("temperature", &ConfigInfo::temperature, "The current simulation temperature.");

	conf_info.def_readonly("current_step", &ConfigInfo::curr_step, "The current time step.");

	conf_info.def_readonly("box", &ConfigInfo::box, "The simulation box, which is an instance of a child class of :class:`BaseBox`.");

}

#endif /* OXPY_BINDINGS_INCLUDES_CONFIGINFO_H_ */
