/*
 * LTCoordination.h
 *
 *  Created on: Oct 17, 2025
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_FORCES_LTCOORDINATION_H_
#define OXPY_BINDINGS_INCLUDES_FORCES_LTCOORDINATION_H_

#include "../../python_defs.h"

#include <Forces/Metadynamics/LTCoordination.h>

void export_LTCoordination(py::module &m) {
	py::class_<LTCoordination, BaseForce, std::shared_ptr<LTCoordination>> force(m, "LTCoordination", R"pbdoc(
A force that depends on a continuous coordination number, as computed from pairs of particles, evaluated from a lookup table.
)pbdoc");

	force.def_readwrite("potential_grid", &LTCoordination::potential_grid, R"c(
The lookup table used to interpolate the potential and the force.
)c");
}

#endif /* OXPY_BINDINGS_INCLUDES_FORCES_LTCOORDINATION_H_ */
