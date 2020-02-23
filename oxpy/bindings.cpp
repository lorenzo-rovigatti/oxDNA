/*
 * bindings.cpp
 *
 *  Created on: Sep 13, 2019
 *      Author: lorenzo
 */

#include "python_defs.h"

#include "OxpyContext.h"
#include "OxpyManager.h"

#include <Particles/BaseParticle.h>
#include <Interactions/BaseInteraction.h>
#include <Utilities/ConfigInfo.h>
#include "vector_matrix_casters.h"

void export_BaseParticle(py::module &m);
void export_IBaseInteraction(py::module &m);
void export_ConfigInfo(py::module &m);

PYBIND11_MODULE(core, m) {
	m.doc() = R"pbdoc(
        Core classes
        ------------

        .. currentmodule:: oxpy.core

        .. autosummary::
           :nosignatures:

           BaseParticle
           ConfigInfo
           IBaseInteraction
           OxpyManager

    )pbdoc";

	export_OxpyContext(m);

	export_SimManager(m);
	export_OxpyManager(m);

	export_BaseParticle(m);
	export_IBaseInteraction(m);
	export_ConfigInfo(m);
}

void export_BaseParticle(py::module &m) {
	pybind11::class_<BaseParticle, std::shared_ptr<BaseParticle>> particle(m, "BaseParticle", R"pbdoc(
        A simulation particle
	)pbdoc");

	particle.def(py::init<>(), R"pbdoc(
		 Default constructor 
	)pbdoc");
	particle.def("is_bonded", &BaseParticle::is_bonded, R"pbdoc(
        Return whether the current particle and arg0 are bonded neighbours

        Parameters
        ----------
        arg0
            The other particle
    )pbdoc");
	particle.def_readwrite("index", &BaseParticle::index, R"pbdoc(
        The index of the particle
    )pbdoc");
	particle.def_readwrite("type", &BaseParticle::type, R"pbdoc(
        The type of the particle
	)pbdoc");
	particle.def_readwrite("btype", &BaseParticle::btype, R"pbdoc(
		The btype of the particle
	)pbdoc");
	particle.def_readwrite("strand_id", &BaseParticle::strand_id, R"pbdoc(
		The id of the strand to which the particle belongs
	)pbdoc");
	particle.def_readwrite("pos", &BaseParticle::pos, R"pbdoc(
		The position of the particle
	)pbdoc");
	particle.def_readwrite("orientation", &BaseParticle::orientation, R"pbdoc(
		The orientation of the particle as a 3x3 matrix
	)pbdoc");
	particle.def_readwrite("vel", &BaseParticle::vel, R"pbdoc(
		The velocity of the particle
	)pbdoc");
	particle.def_readwrite("L", &BaseParticle::L, R"pbdoc(
		The angular momentum of the particle
	)pbdoc");
	particle.def_readwrite("force", &BaseParticle::force, R"pbdoc(
		The force exerted on the particle
	)pbdoc");
	particle.def_readwrite("torque", &BaseParticle::torque, R"pbdoc(
		The torque exerted on the particle
	)pbdoc");
	particle.def_readwrite("ext_potential", &BaseParticle::ext_potential, R"pbdoc(
		The potential energy due to the external forces acting on the particle
	)pbdoc");
	particle.def_readwrite("n3", &BaseParticle::n3, pybind11::return_value_policy::reference, R"pbdoc(
		The n3 neighbour
	)pbdoc");
	particle.def_readwrite("n5", &BaseParticle::n5, pybind11::return_value_policy::reference, R"pbdoc(
		The n5 neighbour
	)pbdoc");
}

void export_ConfigInfo(py::module &m) {
	pybind11::class_<ConfigInfo, std::shared_ptr<ConfigInfo>> conf_info(m, "ConfigInfo", R"pbdoc(
		 This singleton object stores all the details of the simulation (particles, neighbour lists, input file, interaction) 
	)pbdoc");

	conf_info.def("N", &ConfigInfo::N, R"pbdoc(
         Return the current number of particles 
	)pbdoc");
	conf_info.def("particles", &ConfigInfo::particles, pybind11::return_value_policy::reference, R"pbdoc(
		 Return a list of all the particles 
	)pbdoc");
	conf_info.def_readwrite("interaction", &ConfigInfo::interaction, R"pbdoc(
		 The simulation's :py:class:`~oxpy.IBaseInteraction` object 
	)pbdoc");
}

// trampoline class for IBaseInteraction
class PyIBaseInteraction : public IBaseInteraction {
public:
	using IBaseInteraction::IBaseInteraction;

	void init() override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				void,
				IBaseInteraction,
				init
		);
	}

	void allocate_particles(std::vector<BaseParticle *> &particles) override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				void,
				IBaseInteraction,
				allocate_particles,
				particles
		);
	}

	void check_input_sanity(std::vector<BaseParticle *> &particles) override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				void,
				IBaseInteraction,
				check_input_sanity,
				particles
		);
	}

	number pair_interaction(BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false) override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				number,
				IBaseInteraction,
				pair_interaction,
				p,
				q,
				r,
				update_forces
		);

		// suppress warnings
		return 0.;
	}

	number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false) override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				number,
				IBaseInteraction,
				pair_interaction_bonded,
				p,
				q,
				r,
				update_forces
		);

		// suppress warnings
		return 0.;
	}

	number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false) override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				number,
				IBaseInteraction,
				pair_interaction_nonbonded,
				p,
				q,
				r,
				update_forces
		);

		// suppress warnings
		return 0.;
	}

	number pair_interaction_term(int name, BaseParticle *p, BaseParticle *q, LR_vector *r = NULL, bool update_forces = false) override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				number,
				IBaseInteraction,
				pair_interaction_term,
				name,
				p,
				q,
				r,
				update_forces
		);

		// suppress warnings
		return 0.;
	}

	std::map<int, number> get_system_energy_split(std::vector<BaseParticle *> &particles, BaseList *lists) override {
		using ret_type = std::map<int, number>;
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				ret_type,
				IBaseInteraction,
				get_system_energy_split,
				particles,
				lists
		);

		// suppress warnings
		return std::map<int, number>();
	}
};

void export_IBaseInteraction(py::module &m) {
	pybind11::class_<IBaseInteraction, PyIBaseInteraction, std::shared_ptr<IBaseInteraction>> interaction(m, "IBaseInteraction", R"pbdoc(
		 The object that takes care of computing the interaction between the particles 
	)pbdoc");

	interaction.def(py::init<>(), R"pbdoc(
		 Default constructor 
	)pbdoc");
}
