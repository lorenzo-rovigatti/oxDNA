/*
 * Molecule.h
 *
 *  Created on: Feb 14, 2021
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_MOLECULE_H_
#define OXPY_BINDINGS_INCLUDES_MOLECULE_H_

#include "../python_defs.h"

#include <Particles/Molecule.h>

void export_Molecule(py::module &m) {
	py::class_<Molecule, std::shared_ptr<Molecule>> molecule(m, "Molecule", R"pbdoc(
		A set of particles connected by bonds. For DNA and RNA this would be a strand.
	)pbdoc");

	molecule.def_readonly("particles", &Molecule::particles, R"pbdoc(
		The list of particles that compose the molecule.
	)pbdoc");

	molecule.def_property_readonly("id", &Molecule::get_id, R"pbdoc(
		The index of the molecule.
	)pbdoc");

	molecule.def_property_readonly("topology_id", &Molecule::get_topology_id, R"pbdoc(
		The "topology" index of the molecule. For DNA and RNA this value is the id of the strand as listed in the topology file.
	)pbdoc");
}

#endif /* OXPY_BINDINGS_INCLUDES_MOLECULE_H_ */
