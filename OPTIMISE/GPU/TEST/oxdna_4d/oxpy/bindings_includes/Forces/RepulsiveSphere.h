/*
 * RepulsiveSphere.h
 *
 *  Created on: Apr 7, 2022
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_FORCES_REPULSIVESPHERE_H_
#define OXPY_BINDINGS_INCLUDES_FORCES_REPULSIVESPHERE_H_

#include "../../python_defs.h"

#include <Forces/RepulsiveSphere.h>

void export_RepulsiveSphere(py::module &m) {
	py::class_<RepulsiveSphere, BaseForce, std::shared_ptr<RepulsiveSphere>> force(m, "RepulsiveSphere", R"pbdoc(
A confining repulsive sphere.
)pbdoc");

	force.def_readwrite("center", &RepulsiveSphere::_center, R"c(
The center of the sphere.
)c");

	force.def_readwrite("r0", &RepulsiveSphere::_r0, R"c(
The initial radius of the sphere.
)c");

	force.def_readwrite("r_ext", &RepulsiveSphere::_r_ext, R"c(
The radius beyond which the force is not felt any more.
)c");

	force.def_readwrite("rate", &RepulsiveSphere::_rate, R"c(
The rate :math:`A` with which the sphere radius grows over time with the relation :math:`r(t) = r_0 + At`, where :math:`t` is the current time step.
)c");
}

#endif /* OXPY_BINDINGS_INCLUDES_FORCES_REPULSIVESPHERE_H_ */
