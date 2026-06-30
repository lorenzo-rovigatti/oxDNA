/*
 * ParticleStress.h
 *
 *  Created on: Jun 30, 2026
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_OBSERVABLES_PARTICLE_STRESS_H_
#define OXPY_BINDINGS_INCLUDES_OBSERVABLES_PARTICLE_STRESS_H_

#include "../../python_defs.h"

#include <Observables/ParticleStress.h>

void export_ParticleStress(py::module &m) {
	py::class_<ParticleStress, BaseObservable, std::shared_ptr<ParticleStress>> obs(m, "ParticleStress", R"pbdoc(
An observable that computes the stress tensor for each particle.
)pbdoc");

	obs.def(py::init<>(), R"pbdoc(
		The default constructor takes no parameters.
	)pbdoc");

	obs.def("stress_tensors", &ParticleStress::stress_tensors, R"c(
The stress tensors for each particle.

Returns
-------
	List(StressTensor)
		The list of stress tensors for each particle
)c");
}

#endif /* OXPY_BINDINGS_INCLUDES_OBSERVABLES_PARTICLE_STRESS_H_ */
