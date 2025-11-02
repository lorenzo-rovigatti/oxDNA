/*
 * LTCOMTrap.h
 *
 *  Created on: Apr 13, 2022
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_FORCES_LTATANCOMTRAP_H_
#define OXPY_BINDINGS_INCLUDES_FORCES_LTATANCOMTRAP_H_

#include "../../python_defs.h"

#include <Forces/Metadynamics/LTAtanCOMTrap.h>

void export_LTAtanCOMTrap(py::module &m) {
	py::class_<LTAtanCOMTrap, BaseForce, std::shared_ptr<LTAtanCOMTrap>> force(m, "LTAtanCOMTrap", R"pbdoc(
A force depending on the tangent of the ratio of the distances between the centres of mass of groups of particles evaluated from a lookup table.
)pbdoc");

	force.def_readwrite("potential_grid", &LTAtanCOMTrap::potential_grid, R"c(
The lookup table used to interpolate the potential and the force.
)c");
}


#endif /* OXPY_BINDINGS_INCLUDES_FORCES_LTCOMTRAP_H_ */
