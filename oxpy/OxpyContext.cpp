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

OxpyContext::OxpyContext(bool print_coda) :
				_print_coda(print_coda) {

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

	if(_print_coda) {
		OX_LOG(Logger::LOG_NOTHING, "");
		OX_LOG(Logger::LOG_NOTHING, R"c(Please cite at least the following publication for any work that uses the oxDNA simulation package
  * E. Poppleton et al., J. Open Source Softw. 8, 4693 (2023), DOI: 10.21105/joss.04693)c");
	}
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

	context.def(py::init<bool>(), py::arg("print_coda") = true, R"pbdoc(
The default constructor takes an optional parameter that controls whether the final message gets printed or not.

Parameters
----------
print_coda : bool
    If True (default value) a message asking to cite the oxDNA JOSS paper will be printed at the end of the simulation.
	)pbdoc");
	context.def("__enter__", &OxpyContext::enter);
	context.def("__exit__", &OxpyContext::exit);
}
