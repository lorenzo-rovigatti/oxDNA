/*
 * bindings.cpp
 *
 *  Created on: Sep 13, 2019
 *      Author: lorenzo
 */

#include "python_defs.h"

#include "OxpyContext.h"
#include "OxpyManager.h"

PYBIND11_MODULE(oxpy, m) {
	export_OxpyContext(m);

	export_SimManager(m);
	export_OxpyManager(m);
}
