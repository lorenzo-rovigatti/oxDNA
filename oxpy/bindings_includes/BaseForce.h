/*
 * BaseForce.h
 *
 *  Created on: Mar 21, 2020
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_BASEFORCE_H_
#define OXPY_BINDINGS_INCLUDES_BASEFORCE_H_

#include "../python_defs.h"

#include <Forces/BaseForce.h>

// trampoline class for BaseForce
class PyBaseForce : public BaseForce {
public:
	using BaseForce::BaseForce;

	std::tuple<std::vector<int>, std::string> init(input_file &inp, BaseBox *box) override {
		using ret_type = std::tuple<std::vector<int>, std::string>;
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				ret_type,
				BaseForce,
				init,
				inp,
				box
		);

		// suppress warnings
		return std::tuple<std::vector<int>, std::string>();
	}

	LR_vector value(llint step, LR_vector &pos) override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				LR_vector,
				BaseForce,
				value,
				step,
				pos
		);

		// suppress warnings
		return LR_vector();
	}

	number potential(llint step, LR_vector &pos) override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				number,
				BaseForce,
				potential,
				step,
				pos
		);

		// suppress warnings
		return 0.;
	}
};


void export_BaseForce(py::module &m) {
	py::class_<BaseForce, PyBaseForce, std::shared_ptr<BaseForce>> force(m, "BaseForce", R"pbdoc(
		The interface class for forces.
	)pbdoc");

	force.def(py::init<>(), R"pbdoc(
        The default constructor takes no parameters.
	)pbdoc");

	force.def("init", &BaseForce::init, py::arg("inp"), py::arg("box"), R"pbdoc(
        Initialises the force.

        Parameters
        ---------- 
        inp: :class:`input_file`
            The input file of the simulation.
        box: :class:`BaseBox`
            The instance of the simulation's box object.
	)pbdoc");

	force.def("get_group_name", &BaseForce::get_group_name, R"pbdoc(
        Returns the name of the group the force belongs to.

        Returns
        -------
        str
            The name of the group.
	)pbdoc");

	force.def("set_group_name", &BaseForce::set_group_name, py::arg("name"), R"pbdoc(
        Sets the name of the group the force belongs to.

        Parameters
        ---------- 
		name: str
			The new group name.
	)pbdoc");

	force.def("value", &BaseForce::value, py::arg("step"), py::arg("pos"), R"pbdoc(
        Returns the force vector that would act at the given position and time step.

        Parameters
        ----------
        step: int
            The current time step.
        pos: numpy.ndarray
            The position where the force should be computed.

        Returns
        -------
        numpy.ndarray
            The computed force vector.
	)pbdoc");

	force.def("potential", &BaseForce::potential, py::arg("step"), py::arg("pos"), R"pbdoc(
        Returns the potential energy due to the force at the given position and time step.

        Parameters
        ----------
        step: int
            The current time step.
        pos: numpy.ndarray
            The position where the force should be computed.

        Returns
        -------
        float
            The computed potential energy.
	)pbdoc");
}

#endif /* OXPY_BINDINGS_INCLUDES_BASEFORCE_H_ */
