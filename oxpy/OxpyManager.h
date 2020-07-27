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
	OxpyManager(std::string input_filename);
	OxpyManager(input_file input);
	virtual ~OxpyManager();

	std::shared_ptr<ConfigInfo> config_info();
	number system_energy();

	void run(llint steps, bool print_output=true);
};

void export_SimManager(py::module &m);
void export_OxpyManager(py::module &m);

#endif /* OXPY_OXPYMANAGER_H_ */
