/*
 * OxpyManager.h
 *
 *  Created on: Sep 13, 2019
 *      Author: lorenzo
 */

#ifndef OXPY_OXPYMANAGER_H_
#define OXPY_OXPYMANAGER_H_

#include <Managers/SimManager.h>
#include "python_defs.h"

class OxpyManager: public SimManager {
public:
	OxpyManager(std::vector<std::string> args);
	virtual ~OxpyManager();
};

void export_SimManager(py::module &m);
void export_OxpyManager(py::module &m);

#endif /* OXPY_OXPYMANAGER_H_ */
