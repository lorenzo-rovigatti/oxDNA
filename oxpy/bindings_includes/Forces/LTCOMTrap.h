/*
 * LTCOMTrap.h
 *
 *  Created on: Apr 13, 2022
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_FORCES_LTCOMTRAP_H_
#define OXPY_BINDINGS_INCLUDES_FORCES_LTCOMTRAP_H_

#include "../../python_defs.h"

#include <Forces/Metadynamics/LTCOMTrap.h>

void export_LTCOMTrap(py::module &m) {
	py::class_<LTCOMTrap, BaseForce, std::shared_ptr<LTCOMTrap>> force(m, "LTCOMTrap", R"pbdoc(
A force acting between the centres of mass of groups of particles evaluated from a lookup table.
)pbdoc");

	force.def_readwrite("potential_grid", &LTCOMTrap::potential_grid, R"c(
The lookup table used to interpolate the potential and the force.
)c");
}


#endif /* OXPY_BINDINGS_INCLUDES_FORCES_LTCOMTRAP_H_ */
