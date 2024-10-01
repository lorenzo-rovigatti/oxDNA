/*
 * OxpyContext.h
 *
 *  Created on: Sep 17, 2019
 *      Author: lorenzo
 */

#ifndef OXPY_OXPYCONTEXT_H_
#define OXPY_OXPYCONTEXT_H_

#include "python_defs.h"

class OxpyContext {
public:
	OxpyContext(bool print_coda=true);
	virtual ~OxpyContext();

	void enter();
	void exit(py::object exc_type, py::object exc_value, py::object traceback);

private:
	bool _print_coda;
};

void export_OxpyContext(py::module &m);

#endif /* OXPY_OXPYCONTEXT_H_ */
