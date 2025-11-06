/*
 * BaseParticle.h
 *
 *  Created on: Mar 21, 2020
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_BASEPARTICLE_H_
#define OXPY_BINDINGS_INCLUDES_BASEPARTICLE_H_

#include "../python_defs.h"

#include <Particles/BaseParticle.h>

void export_BaseParticle(py::module &m) {
	py::class_<BaseParticle, std::shared_ptr<BaseParticle>> particle(m, "BaseParticle", R"pbdoc(
        A simulation particle.
	)pbdoc");

	particle.def(py::init<>(), R"pbdoc(
		 The default constructor takes no parameters.
	)pbdoc");
	particle.def("is_bonded", &BaseParticle::is_bonded, py::arg("q"), R"pbdoc(
        Return whether the current particle and q are bonded neighbours.

        Parameters
        ----------
        q: :class:`BaseParticle`
            The other Particle.

        Returns
        -------
        bool
            True if the current particle and :attr:`q` are bonded neighbours.
    )pbdoc");
	particle.def_readwrite("index", &BaseParticle::index, R"pbdoc(
        The index of the particle.
    )pbdoc");
	particle.def_readwrite("type", &BaseParticle::type, R"pbdoc(
        The type of the particle.
	)pbdoc");
	particle.def_readwrite("btype", &BaseParticle::btype, R"pbdoc(
		The btype of the particle.
	)pbdoc");
	particle.def_readwrite("strand_id", &BaseParticle::strand_id, R"pbdoc(
		The id of the strand to which the particle belongs.
	)pbdoc");
	particle.def_readwrite("pos", &BaseParticle::pos, R"pbdoc(
		The position of the particle.
	)pbdoc");
	particle.def_readwrite("orientation", &BaseParticle::orientation, R"pbdoc(
		The orientation of the particle as a 3x3 matrix.
	)pbdoc");
	particle.def_readwrite("vel", &BaseParticle::vel, R"pbdoc(
		The velocity of the particle.
	)pbdoc");
	particle.def_readwrite("L", &BaseParticle::L, R"pbdoc(
		The angular momentum of the particle.
	)pbdoc");
	particle.def_readwrite("force", &BaseParticle::force, R"pbdoc(
		The force exerted on the particle.
	)pbdoc");
	particle.def_readwrite("torque", &BaseParticle::torque, R"pbdoc(
		The torque exerted on the particle.
	)pbdoc");
	particle.def_readwrite("ext_potential", &BaseParticle::ext_potential, R"pbdoc(
		The potential energy due to the external forces acting on the particle.
	)pbdoc");
	particle.def_readwrite("n3", &BaseParticle::n3, py::return_value_policy::reference, R"pbdoc(
		The n3 neighbour.
	)pbdoc");
	particle.def_readwrite("n5", &BaseParticle::n5, py::return_value_policy::reference, R"pbdoc(
		The n5 neighbour.
	)pbdoc");
	particle.def_readwrite("int_centers", &BaseParticle::int_centers, R"pbdoc(
		The interaction sites of the particle.
	)pbdoc");
}

#endif /* OXPY_BINDINGS_INCLUDES_BASEPARTICLE_H_ */
