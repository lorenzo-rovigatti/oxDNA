/*
 * MovingTrap.h
 *
 *  Created on: Apr 20, 2022
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_FORCES_MOVINGTRAP_H_
#define OXPY_BINDINGS_INCLUDES_FORCES_MOVINGTRAP_H_

#include "../../python_defs.h"

#include <Forces/MovingTrap.h>

void export_MovingTrap(py::module &m) {
	py::class_<MovingTrap, BaseForce, std::shared_ptr<MovingTrap>> force(m, "MovingTrap", R"pbdoc(
A harmonic force computed by using as a reference a point that optionally move with a constant velocity along a direction.
)pbdoc");

	force.def_readwrite("dir", &MovingTrap::_direction, R"c(
The direction towards which the reference point moves.
)c");

	force.def_readwrite("pos0", &MovingTrap::_pos0, R"c(
The initial position of the reference point.
)c");

	force.def_readwrite("stiff", &MovingTrap::_stiff, R"c(
The stiffness of the harmonic force.
)c");

	force.def_readwrite("rate", &MovingTrap::_rate, R"c(
The velocity (in units of time steps) with which the reference point moves. If set to 0 the point will remain in its initial position.
)c");
}

#endif /* OXPY_BINDINGS_INCLUDES_FORCES_MOVINGTRAP_H_ */
