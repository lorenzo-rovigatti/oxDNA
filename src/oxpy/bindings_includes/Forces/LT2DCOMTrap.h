/*
 * LT2DCOMTrap.h
 *
 *  Created on: Apr 13, 2022
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_FORCES_LT2DCOMTRAP_H_
#define OXPY_BINDINGS_INCLUDES_FORCES_LT2DCOMTRAP_H_

#include "../../python_defs.h"

#include <Forces/Metadynamics/LT2DCOMTrap.h>

void export_LT2DCOMTrap(py::module &m) {
	py::class_<LT2DCOMTrap, BaseForce, std::shared_ptr<LT2DCOMTrap>> force(m, "LT2DCOMTrap", R"pbdoc(
A force acting between the centres of mass of groups of particles evaluated from a 2D lookup table.

The force is evaluated by using two collective variables that are the distances between two pairs of centres of mass. 
)pbdoc");

	force.def_readwrite("potential_grid", &LT2DCOMTrap::potential_grid, R"c(
The lookup table used to interpolate the potential and the force.
)c");
}

#endif /* OXPY_BINDINGS_INCLUDES_FORCES_LT2DCOMTRAP_H_ */
