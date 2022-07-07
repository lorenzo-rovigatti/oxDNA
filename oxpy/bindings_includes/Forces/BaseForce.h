/*
 * BaseForce.h
 *
 *  Created on: Mar 21, 2020
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_FORCES_BASEFORCE_H_
#define OXPY_BINDINGS_INCLUDES_FORCES_BASEFORCE_H_

#include "../../python_defs.h"

#include "LTCOMTrap.h"
#include "LTCOMAngleTrap.h"
#include "LTAtanCOMTrap.h"
#include "LT2DCOMTrap.h"
#include "MovingTrap.h"
#include "RepulsiveSphere.h"

#include <Forces/BaseForce.h>

// trampoline class for BaseForce
class PyBaseForce : public BaseForce {
public:
	using BaseForce::BaseForce;

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override {
		using ret_type = std::tuple<std::vector<int>, std::string>;
		PYBIND11_OVERLOAD( // @suppress("Unused return value")
				ret_type,
				BaseForce,
				init,
				inp
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
	py::module sub_m = m.def_submodule("forces");

	py::class_<BaseForce, PyBaseForce, std::shared_ptr<BaseForce>> force(sub_m, "BaseForce", R"pbdoc(
		The interface class for forces.
	)pbdoc");

	force.def(py::init<>(), R"pbdoc(
        The default constructor takes no parameters.
	)pbdoc");

	force.def("init", &BaseForce::init, py::arg("inp"), R"pbdoc(
        Initialise the force.

        Parameters
        ---------- 
        inp: :class:`input_file`
            The input file of the simulation.
	)pbdoc");

	force.def_property("group_name", &BaseForce::get_group_name, &BaseForce::set_group_name, "The name of the group the force belongs to.");
	force.def_property("id", &BaseForce::get_id, &BaseForce::set_id, "The id of the force.");
	force.def_property_readonly("type", &BaseForce::get_type, "the string used in the external force file to specify the force type (*e.g.* `trap`).");

	force.def("value", &BaseForce::value, py::arg("step"), py::arg("pos"), R"pbdoc(
        Return the force vector that would act at the given position and time step.

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
        Return the potential energy due to the force at the given position and time step.

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

	force.def_readwrite("stiff", &BaseForce::_stiff, R"pbdoc(
The stiffness (= strength) of the force. The particular meaning of this attribute depends on the subclass external force.
)pbdoc");

	force.def_readwrite("rate", &BaseForce::_rate, R"pbdoc(
The rate of change of the force. The particular meaning of this attribute depends on the subclass external force.
)pbdoc");

	force.def_readwrite("pos0", &BaseForce::_pos0, R"pbdoc(
The reference position associated to the force. The particular meaning of this attribute depends on the subclass external force.
)pbdoc");

	force.def("as_RepulsiveSphere", [](BaseForce &f){ return dynamic_cast<RepulsiveSphere *>(&f); }, R"pbdoc(
Attempt to cast the current force as a :py:class:`RepulsiveSphere` object and return it.

Returns
-------
	:py:class:`RepulsiveSphere`
		The force cast as a :py:class:`RepulsiveSphere` or `None` if the casting fails.
)pbdoc");

	export_LTCOMTrap(sub_m);
	export_LTCOMAngleTrap(sub_m);
	export_LTAtanCOMTrap(sub_m);
	export_LT2DCOMTrap(sub_m);
	export_MovingTrap(sub_m);
	export_RepulsiveSphere(sub_m);
}

#endif /* OXPY_BINDINGS_INCLUDES_FORCES_BASEFORCE_H_ */
