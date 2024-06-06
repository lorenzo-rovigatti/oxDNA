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

Options can also be deleted::

	del my_input["log_file"]

To check if a specific option has been set you can use `in`::

	if "debug" in my_input:
		# do something
	)pbdoc");

	input.def(py::init<>());

	input.def("__getitem__", &input_file::get_value,
		py::arg("key") = '\0',
		py::arg("mandatory") = 0,
		py::arg("found") = true
	);

	input.def("__setitem__", &input_file::set_value);

	input.def("__delitem__", &input_file::unset_value);

	input.def("__str__", &input_file::to_string);

	// this enables the use of "if 'something' in input_file"
	input.def("__contains__", [](input_file &inp, std::string v) {
		try {
			bool found;
			inp.get_value(v, 1, found);
			return true;
		}
		catch (const oxDNAException &e){
			return false;
		}
	});

	input.def("init_from_filename", &input_file::init_from_filename, py::arg("filename"), R"pbdoc(
        Initialise the object from the input file passed as parameter.

        Parameters
        ----------
        filename: str
            The name of the input file whence the options will be loaded.
	)pbdoc");

	input.def("get_bool", [](input_file &inp, std::string &key) {
		bool found;
		if (inp.true_values.find(inp.get_value(key, 0, found)) != inp.true_values.end()) {
			return true;
		}
		else {
			if (inp.false_values.find(inp.get_value(key, 0, found)) != inp.false_values.end()){
				return false;
			}
			else {
				throw std::invalid_argument("boolean key " + key + " is invalid");
			}
		}
	}, R"pbdoc(
		Return the boolean value of an input option.

		Returns
		-------
		bool
			The boolean value of the option.
	)pbdoc");

	input.def("keys", [](input_file &inp) {
		std::vector<std::string> keys;
		for (auto &it : inp.keys) {
			keys.push_back(it.first);
		}
		return keys;
	}, R"pbdoc(
		Return the list of keys stored in the input file. For instance, the following code prints all the options stored in an input file::

			keys = my_input.keys()
			for key in keys:
				print(f"Option '{key}' is set to '{my_input[key]}'")

		Returns
		-------
		List(str)
			The list of keys stored in the input file.
	)pbdoc");
}

#endif /* OXPY_BINDINGS_INCLUDES_INPUT_FILE_H_ */
