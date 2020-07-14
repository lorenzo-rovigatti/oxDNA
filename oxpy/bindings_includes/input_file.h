/*
 * input_file.h
 *
 *  Created on: Jun 20, 2020
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_INPUT_FILE_H_
#define OXPY_BINDINGS_INCLUDES_INPUT_FILE_H_

#include "../python_defs.h"

#include <Utilities/parse_input/parse_input.h>

void export_input_file(py::module &m) {
	py::class_<input_file, std::shared_ptr<input_file>> input(m, "InputFile", R"pbdoc(
        Handles the options provided by the user to control and tune the behaviour of the simulation/analysis/generation of the initial configuration.

        Options can be read, set or overwritten by using square brackets::

            my_input["sim_type"] = "VMMC"
            print(my_input["sim_type"])
	)pbdoc");

	input.def(py::init<>());
	input.def("__getitem__", &input_file::get_value);
	input.def("__setitem__", &input_file::set_value);
	input.def("__str__", &input_file::to_string);
	input.def("init_from_filename", &input_file::init_from_filename, py::arg("filename"), R"pbdoc(
        Initialise the object from the input file passed as parameter.

        Parameters
        ----------
        filename: str
            The name of the input file whence the options will be loaded.
	)pbdoc");
}

#endif /* OXPY_BINDINGS_INCLUDES_INPUT_FILE_H_ */
