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
		OX_LOG(Logger::LOG_NOTHING, R"c(Please cite these publications for any work that uses the oxDNA simulation package
		- for the code:
			* P. Šulc et al., J. Chem. Phys. 137, 135101 (2012)
			* L. Rovigatti et al., J. Comput. Chem. 36, 1 (2015)
		- for the oxDNA model:
			* T. E. Ouldridge et al., J. Chem. Phys, 134, 085101 (2011)
		- for the oxDNA2 model:
			* B. E. K. Snodin et al., J. Chem. Phys. 142, 234901 (2015)
		- for the oxRNA model:
			* P. Šulc et al., J. Chem. Phys. 140, 235102 (2014))c");
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

	context.def(py::init<bool>(), py::arg("print_coda") = true);
	context.def("__enter__", &OxpyContext::enter);
	context.def("__exit__", &OxpyContext::exit);
}
