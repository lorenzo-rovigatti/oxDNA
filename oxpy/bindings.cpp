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

void export_BaseParticle(py::module &m);
void export_IBaseInteraction(py::module &m);
void export_ConfigInfo(py::module &m);

PYBIND11_MODULE(oxpy, m) {
	export_OxpyContext(m);

	export_SimManager(m);
	export_OxpyManager(m);

	export_BaseParticle(m);
	export_IBaseInteraction(m);
	export_ConfigInfo(m);
}

void export_BaseParticle(py::module &m) {
	pybind11::class_<BaseParticle, std::shared_ptr<BaseParticle>> particle(m, "BaseParticle");

	particle
		.def(py::init<>());
}

void export_ConfigInfo(py::module &m) {
	pybind11::class_<ConfigInfo, std::shared_ptr<ConfigInfo>> conf_info(m, "ConfigInfo");

	conf_info
		.def("N", &ConfigInfo::N)
		.def("particles", &ConfigInfo::particles)
		.def_readwrite("interaction", &ConfigInfo::interaction);
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
	pybind11::class_<IBaseInteraction, PyIBaseInteraction, std::shared_ptr<IBaseInteraction>> interaction(m, "IBaseInteraction");

	interaction
		.def(py::init<>());
}
