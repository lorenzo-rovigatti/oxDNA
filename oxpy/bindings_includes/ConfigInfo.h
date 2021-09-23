/*
 * ConfigInfo.h
 *
 *  Created on: Mar 21, 2020
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_CONFIGINFO_H_
#define OXPY_BINDINGS_INCLUDES_CONFIGINFO_H_

#include "../python_defs.h"

#include <Utilities/ConfigInfo.h>

void export_ConfigInfo(py::module &m) {
	py::class_<ConfigInfo, std::shared_ptr<ConfigInfo>> conf_info(m, "ConfigInfo", R"pbdoc(
		 This singleton object stores all the details of the simulation (particles, neighbour lists, input file, interaction) 
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
}

#endif /* OXPY_BINDINGS_INCLUDES_CONFIGINFO_H_ */
