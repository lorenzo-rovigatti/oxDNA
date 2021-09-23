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
}

#endif /* OXPY_BINDINGS_INCLUDES_MOLECULE_H_ */
