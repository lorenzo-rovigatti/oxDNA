/*
 * HBList.h
 *
 *  Created on: Apr 28, 2022
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_OBSERVABLES_HBLIST_H_
#define OXPY_BINDINGS_INCLUDES_OBSERVABLES_HBLIST_H_

#include "../../python_defs.h"

#include <Observables/HBList.h>

void export_HBList(py::module &m) {
	py::class_<HBList, BaseObservable, std::shared_ptr<HBList>> obs(m, "HBList", R"pbdoc(
An observable that computes a list of hydrogen-bonded nucleotides.
)pbdoc");

	obs.def(py::init<>(), R"pbdoc(
		The default constructor takes no parameters.
	)pbdoc");

	obs.def("hb_list", &HBList::hb_list, R"c(
The list of hydrogen-bonded pairs.

Returns
-------
	List(Tuple(int))
		The list of pairs involved in hydrogen bonds
)c");
}

#endif /* OXPY_BINDINGS_INCLUDES_OBSERVABLES_HBLIST_H_ */
