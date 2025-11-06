/*
 * RNANucleotide.h
 *
 *  Created on: Mar 22, 2020
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_RNANUCLEOTIDE_H_
#define OXPY_BINDINGS_INCLUDES_RNANUCLEOTIDE_H_

#include "../python_defs.h"

#include <Particles/RNANucleotide.h>

void export_RNANucleotide(py::module &m) {
	py::class_<RNANucleotide, BaseParticle, std::shared_ptr<RNANucleotide>> nucleotide(m, "RNANucleotide", R"pbdoc(
        A RNA nucleotide. This class exposes the :class:`Site` enum to access the nucleotide's sites.
    )pbdoc");

	nucleotide.def(py::init<>(), R"pbdoc(
        The constructor takes no parameters.
    )pbdoc");

	nucleotide.def("backbone_site", [](RNANucleotide &nucl) {
		return nucl.pos + nucl.int_centers[nucl.BACK];
	}, R"pbdoc(
		Return the position of the backbone site.

		Returns
		-------
		numpy.ndarray
			The position of the site.
	)pbdoc");

	nucleotide.def("stacking_site", [](RNANucleotide &nucl) {
		return nucl.pos + nucl.int_centers[nucl.STACK];
	}, R"pbdoc(
        Returns the position of the stacking site.

        Returns
        -------
        numpy.ndarray
            The position of the site.
	)pbdoc");

	nucleotide.def("base_site", [](RNANucleotide &nucl) {
		return nucl.pos + nucl.int_centers[nucl.BASE];
	}, R"pbdoc(
        Return the position of the base site.

        Returns
        -------
        numpy.ndarray
            The position of the site.
	)pbdoc");

	nucleotide.def("stack3_site", [](RNANucleotide &nucl) {
		return nucl.pos + nucl.int_centers[nucl.STACK_3];
	}, R"pbdoc(
        Return the position of the stack3 site.

        Returns
        -------
        numpy.ndarray
            The position of the site.
	)pbdoc");

	nucleotide.def("stack5_site", [](RNANucleotide &nucl) {
		return nucl.pos + nucl.int_centers[nucl.STACK_5];
	}, R"pbdoc(
        Return the position of the stack5 site.

        Returns
        -------
        numpy.ndarray
            The position of the site.
	)pbdoc");

	nucleotide.def("bbvector3_site", [](RNANucleotide &nucl) {
		return nucl.pos + nucl.int_centers[nucl.BBVECTOR_3];
	}, R"pbdoc(
        Return the position of the bbvector3 site.

        Returns
        -------
        numpy.ndarray
            The position of the site.
	)pbdoc");

	nucleotide.def("bbvector5_site", [](RNANucleotide &nucl) {
		return nucl.pos + nucl.int_centers[nucl.BBVECTOR_5];
	}, R"pbdoc(
        Return the position of the bbvector5 site.

        Returns
        -------
        numpy.ndarray
            The position of the site.
	)pbdoc");
}

#endif /* OXPY_BINDINGS_INCLUDES_BASEPARTICLE_H_ */
