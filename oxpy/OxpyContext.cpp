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
#include <Utilities/ConfigInfo.h>

OxpyContext::OxpyContext() {

}

OxpyContext::~OxpyContext() {

}

void OxpyContext::enter() {
	try {
		Logger::init();
	}
	catch (oxDNAException &e) {
		// the Logger object may be initialised only once
	}
	SignalManager::manage_segfault();
	TimingManager::init();
}

void OxpyContext::exit(py::object exc_type, py::object exc_value, py::object traceback) {
	TimingManager::clear();
	ConfigInfo::clear();
}

void export_OxpyContext(py::module &m) {
	pybind11::class_<OxpyContext, std::shared_ptr<OxpyContext>> context(m, "Context", R"pbdoc(
        A `context manager <https://book.pythontips.com/en/latest/context_managers.html>`_ that 
        performs some initialisation and clean up. 

        Use it with a :code:`with` command before using
        any of the :code:`oxpy` functionalities::

            with oxpy.Context():
                # your code using oxpy goes here
	)pbdoc");

	context
		.def(py::init<>())
		.def("__enter__", &OxpyContext::enter)
		.def("__exit__", &OxpyContext::exit);
}
