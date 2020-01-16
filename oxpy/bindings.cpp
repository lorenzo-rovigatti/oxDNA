/*
 * bindings.cpp
 *
 *  Created on: Sep 13, 2019
 *      Author: lorenzo
 */

#include "python_defs.h"

#include "OxpyContext.h"
#include "OxpyManager.h"

#include <Particles/BaseParticle.h>
#include <Interactions/BaseInteraction.h>
#include <Utilities/ConfigInfo.h>

void export_BaseParticle(py::module &m);
void export_IBaseInteraction(py::module &m);
void export_ConfigInfo(py::module &m);

PYBIND11_MODULE(oxpy, m) {
	export_OxpyContext(m);

	export_SimManager(m);
	export_OxpyManager(m);

	export_BaseParticle(m);
//	export_IBaseInteraction(m);
	export_ConfigInfo(m);
}

void export_BaseParticle(py::module &m) {
	pybind11::class_<BaseParticle, std::shared_ptr<BaseParticle>> particle(m, "BaseParticle");

	particle
		.def(py::init<>());
}

void export_IBaseInteraction(py::module &m) {
	pybind11::class_<IBaseInteraction, std::shared_ptr<IBaseInteraction>> interaction(m, "IBaseInteraction");

	interaction
		.def(py::init<>());
}

void export_ConfigInfo(py::module &m) {
	pybind11::class_<ConfigInfo, std::shared_ptr<ConfigInfo>> conf_info(m, "ConfigInfo");

	conf_info
//		.def_readwrite("particles", &ConfigInfo::particles)
		.def_readwrite("interaction", &ConfigInfo::interaction)
		.def_readwrite("N", &ConfigInfo::N);
}
