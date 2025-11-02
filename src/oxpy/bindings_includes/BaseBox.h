/*
 * BaseBox.h
 *
 *  Created on: Jun 14, 2022
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_FORCES_BASEBOX_H_
#define OXPY_BINDINGS_INCLUDES_FORCES_BASEBOX_H_

#include "../python_defs.h"

#include <Boxes/BaseBox.h>

// trampoline class for BaseBox
class PyBaseBox : public BaseBox {
public:
	using BaseBox::BaseBox;

	void init(number Lx, number Ly, number Lz) override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				void,
				BaseBox,
				init,
				Lx,
				Ly,
				Lz
		);
	}

	void get_settings(input_file &inp) override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				void,
				BaseBox,
				get_settings,
				inp
		);
	}

	LR_vector min_image(const LR_vector &v1, const LR_vector &v2) const override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				LR_vector,
				BaseBox,
				min_image,
				v1,
				v2
		);

		return LR_vector();
	}

	number sqr_min_image_distance(const LR_vector &v1, const LR_vector &v2) const override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				number,
				BaseBox,
				sqr_min_image_distance,
				v1,
				v2
		);

		return 0.;
	}

	LR_vector normalised_in_box(const LR_vector &v) override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				LR_vector,
				BaseBox,
				normalised_in_box,
				v
		);

		return LR_vector();
	}

	LR_vector box_sides() const override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				LR_vector &,
				BaseBox,
				box_sides
		);

		return LR_vector();
	}

	number V() override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				number,
				BaseBox,
				V
		);

		return 0.;
	}

	LR_vector get_abs_pos(BaseParticle *p) override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				LR_vector,
				BaseBox,
				get_abs_pos
				p
		);

		return LR_vector();
	}

	void shift_particle(BaseParticle *p, LR_vector &amount) override {
		PYBIND11_OVERLOAD_PURE( // @suppress("Unused return value")
				void,
				BaseBox,
				shift_particle,
				p,
				amount
		)
	}
};

void export_BaseBox(py::module &m) {
	py::class_<BaseBox, PyBaseBox, std::shared_ptr<BaseBox>> box(m, "BaseBox", R"pbdoc(
		The interface class for forces.
	)pbdoc");

	box.def(py::init<>(), R"pbdoc(
        The default constructor takes no parameters.
	)pbdoc");

	box.def("get_settings", &BaseBox::get_settings, py::arg("inp"), R"pbdoc(
        Parse the options relative to this box.

        Parameters
        ---------- 
        inp: :class:`input_file`
            The input file of the simulation.
	)pbdoc");

	box.def("init", &BaseBox::init, R"pbdoc(
        Initialise the box. Calling this method fires off a "box_changed" event.
	)pbdoc");

	box.def_property_readonly("box_sides", &BaseBox::box_sides);

	box.def("min_image", static_cast<LR_vector (BaseBox::*)(const BaseParticle &, const BaseParticle &)>(&BaseBox::min_image), py::arg("p"), py::arg("q"), R"pbdoc(
		Compute and return the PBC-aware distance vector between two particles. The resulting vector points from particle `p` to particle `q`.

		Parameters
		----------
		p: :class:`BaseParticle`
			The first particle.
		q: :class:`BaseParticle`
			The second particle.
		Returns
		-------
		numpy.ndarray
			The distance vector connecting `p` to `q`.
	)pbdoc");

	box.def("sqr_min_image_distance", static_cast<number (BaseBox::*)(const BaseParticle &, const BaseParticle &)>(&BaseBox::sqr_min_image_distance), py::arg("p"), py::arg("q"), R"pbdoc(
		Compute and return the PBC-aware squared distance between two particles. The result is the squared length of the vector that points from particle `p` to particle `q`.

		Parameters
		----------
		p: :class:`BaseParticle`
			The first particle.
		q: :class:`BaseParticle`
			The second particle.
		Returns
		-------
		float
			The squared length of the vector connecting `p` to `q`. 
		
	)pbdoc");

	box.def_property_readonly("V", &BaseBox::V, R"pbdoc(
		The volume of the simulation box.
	)pbdoc");

	box.def("get_abs_pos", &BaseBox::get_abs_pos, py::arg("p"), R"pbdoc(
		Return the absolute (unwrapped) position of a particle.

		Parameters
		----------
		p: :class:`BaseParticle`
			The particle of interest.
		Returns
		-------
		numpy.ndarray
			The absolute (unwrapped) position of the given particle.
	)pbdoc");

}

#endif /* OXPY_BINDINGS_INCLUDES_FORCES_BASEBOX_H_ */
