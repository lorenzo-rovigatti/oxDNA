/*
 * BaseInteraction.h
 *
 *  Created on: Mar 21, 2020
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_BASEINTERACTION_H_
#define OXPY_BINDINGS_INCLUDES_BASEINTERACTION_H_

#include "../python_defs.h"

#include <Interactions/BaseInteraction.h>

// trampoline class for BaseInteraction
class PyBaseInteraction : public BaseInteraction {
public:
	using BaseInteraction::BaseInteraction;

	void init() override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				void,
				BaseInteraction,
				init
		);
	}

	void allocate_particles(std::vector<BaseParticle *> &particles) override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				void,
				BaseInteraction,
				allocate_particles,
				particles
		);
	}

	void check_input_sanity(std::vector<BaseParticle *> &particles) override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				void,
				BaseInteraction,
				check_input_sanity,
				particles
		);
	}

	number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false) override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				number,
				BaseInteraction,
				pair_interaction,
				p,
				q,
				compute_r,
				update_forces
		);

		// suppress warnings
		return 0.;
	}

	number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false) override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				number,
				BaseInteraction,
				pair_interaction_bonded,
				p,
				q,
				compute_r,
				update_forces
		);

		// suppress warnings
		return 0.;
	}

	number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false) override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				number,
				BaseInteraction,
				pair_interaction_nonbonded,
				p,
				q,
				compute_r,
				update_forces
		);

		// suppress warnings
		return 0.;
	}

	number pair_interaction_term(int name, BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false) override {
		PYBIND11_OVERLOAD( // @suppress("Unused return value")
				number,
				BaseInteraction,
				pair_interaction_term,
				name,
				p,
				q,
				compute_r,
				update_forces
		);

		// suppress warnings
		return 0.;
	}

	std::map<int, number> get_system_energy_split(std::vector<BaseParticle *> &particles, BaseList *lists) override {
		using ret_type = std::map<int, number>;
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				ret_type,
				BaseInteraction,
				get_system_energy_split,
				particles,
				lists
		);

		// suppress warnings
		return std::map<int, number>();
	}
};

void export_BaseInteraction(py::module &m) {
	py::class_<BaseInteraction, PyBaseInteraction, std::shared_ptr<BaseInteraction>> interaction(m, "BaseInteraction", R"pbdoc(
		The class that takes care of computing the interaction between the particles.
	)pbdoc");

	interaction.def(py::init<>(), R"pbdoc(
		The default constructor takes no parameters. 
	)pbdoc");

	interaction.def("set_computed_r", &BaseInteraction::set_computed_r, py::arg("r"), R"pbdoc(
        Set the distance vector used by the `pair_interaction_*` methods when they are called with `compute_r = False` (see :meth:`pair_interaction` for additional details).

        Parameters
        ----------
        r : numpy.ndarray
            The distance vector to be stored.
	)pbdoc");

	interaction.def("pair_interaction", &BaseInteraction::pair_interaction, py::arg("p"), py::arg("q"), py::arg("compute_r") = true, py::arg("update_forces") = false, R"pbdoc(
        Compute the pair interaction between p and q.

        Parameters
        ----------
        p : :class:`BaseParticle`
            The first particle of the pair. Note that some interactions require that the two particles are passed to the method with a specific order.
        q : :class:`BaseParticle`
            The second particle of the pair.
        compute_r : bool
            If True (default value), the distance between :attr:`p` and :attr:`q` will be computed from scratch. If not, it will use a private member that can be
            set through the :meth:`set_computed_r` method.
        update_forces : bool
            If True, the forces and torques acting on the two particles will be updated (defaults to False).

        Returns
        -------
        float
            The energy of the pair interaction.
    )pbdoc");

	interaction.def("pair_interaction_bonded", &BaseInteraction::pair_interaction_bonded, py::arg("p"), py::arg("q"), py::arg("compute_r") = true, py::arg("update_forces") = false, R"pbdoc(
        Compute the bonded pair interaction between p and q. See :meth:`pair_interaction` for details on the parameters and on the return value.
	)pbdoc");

	interaction.def("pair_interaction_nonbonded", &BaseInteraction::pair_interaction_nonbonded, py::arg("p"), py::arg("q"), py::arg("compute_r") = true, py::arg("update_forces") = false, R"pbdoc(
        Compute the unbonded pair interaction between p and q. See :meth:`pair_interaction` for details on the parameters and on the return value.
	)pbdoc");

	interaction.def("pair_interaction_term", &BaseInteraction::pair_interaction_term, py::arg("term_id"), py::arg("p"), py::arg("q"), py::arg("compute_r") = true, py::arg("update_forces") = false, R"pbdoc(
		Compute the pair interaction between p and q due to the specified term.

		Parameters
		----------
		term_id : int
			The id of the energy term to be computed.
		p : :class:`BaseParticle`
			The first particle of the pair. Note that some interactions require that the two particles are passed to the method with a specific order.
		q : :class:`BaseParticle`
			The second particle of the pair.
		compute_r : bool
			If True (default value), the distance between :attr:`p` and :attr:`q` will be computed from scratch. If not, it will use a private member that can be
			set through the :meth:`set_computed_r` method.
		update_forces : bool
			If True, the forces and torques acting on the two particles will be updated (defaults to False).

		Returns
		-------
		float
			The energy of the pair interaction due to the specified term.
	)pbdoc");

	interaction.def("has_custom_stress_tensor", &BaseInteraction::has_custom_stress_tensor, R"pbdoc(
		Return whether the interaction computes the stress tensor internally or not.

		Returns
		-------
		bool
			True if the interaction computes the stress tensor internally, False otherwise.
    )pbdoc");

	interaction.def("begin_energy_computation", &BaseInteraction::begin_energy_computation, R"pbdoc(
		Signals the interaction that an energy (or force) computation is about to begin.
    )pbdoc");
}

#endif /* OXPY_BINDINGS_INCLUDES_BASEINTERACTION_H_ */
