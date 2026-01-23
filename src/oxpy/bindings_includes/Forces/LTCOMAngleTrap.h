/*
 * LTCOMAngleTrap.h
 *
 *  Created on: May 06, 2022
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_FORCES_LTCOMANGLETRAP_H_
#define OXPY_BINDINGS_INCLUDES_FORCES_LTCOMANGLETRAP_H_

#include "../../python_defs.h"

#include <Forces/Metadynamics/LTCOMAngleTrap.h>

void export_LTCOMAngleTrap(py::module &m) {
	py::class_<LTCOMAngleTrap, BaseForce, std::shared_ptr<LTCOMAngleTrap>> force(m, "LTCOMAngleTrap", R"pbdoc(
A force that depends on the angle between the distances connecting the centres of mass of three sets of nucleotides.
)pbdoc");

	force.def_readwrite("potential_grid", &LTCOMAngleTrap::potential_grid, R"c(
The lookup table used to interpolate the potential and the force.
)c");
}


#endif /* OXPY_BINDINGS_INCLUDES_FORCES_LTCOMTRAP_H_ */
