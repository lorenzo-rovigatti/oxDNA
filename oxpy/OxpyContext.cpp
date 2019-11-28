/*
 * OxpyContext.cpp
 *
 *  Created on: Sep 17, 2019
 *      Author: lorenzo
 */

#include "OxpyContext.h"

#include <Utilities/Logger.h>
#include <Utilities/SignalManager.h>
#include <Utilities/Timings.h>

OxpyContext::OxpyContext() {

}

OxpyContext::~OxpyContext() {

}

void OxpyContext::enter() {
	Logger::init();
	SignalManager::manage_segfault();
	TimingManager::init();
}

void OxpyContext::exit(py::object exc_type, py::object exc_value, py::object traceback) {
	TimingManager::clear();
}

void export_OxpyContext(py::module &m) {
	pybind11::class_<OxpyContext, std::shared_ptr<OxpyContext>> context(m, "Context");

	context
		.def(py::init<>())
		.def("__enter__", &OxpyContext::enter)
		.def("__exit__", &OxpyContext::exit);
}
