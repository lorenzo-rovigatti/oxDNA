/*
 * DNANucleotide.h
 *
 *  Created on: Mar 22, 2020
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_DNANUCLEOTIDE_H_
#define OXPY_BINDINGS_INCLUDES_DNANUCLEOTIDE_H_

#include "../python_defs.h"

#include <Particles/DNANucleotide.h>

void export_DNANucleotide(py::module &m) {
	py::class_<DNANucleotide, BaseParticle, std::shared_ptr<DNANucleotide>> nucleotide(m, "DNANucleotide", R"pbdoc(
        A DNA nucleotide.
    )pbdoc");

	nucleotide.def(py::init<bool>(), py::arg("grooving"), R"pbdoc(
        The constructor takes a single mandatory parameter.

        Parameters
        ----------
        grooving: bool
            The type of DNA to consider (with major/minor grooving or without).
    )pbdoc");

	nucleotide.def("backbone_site", [](DNANucleotide &nucl) {
		return nucl.pos + nucl.int_centers[nucl.BACK];
	}, R"pbdoc(
        Return the position of the backbone site.

        Returns
        -------
        numpy.ndarray
            The position of the site.
	)pbdoc");

	nucleotide.def("stacking_site", [](DNANucleotide &nucl) {
		return nucl.pos + nucl.int_centers[nucl.STACK];
	}, R"pbdoc(
    Return the position of the stacking site.

    Returns
    -------
    numpy.ndarray
        The position of the site.
	)pbdoc");

	nucleotide.def("base_site", [](DNANucleotide &nucl) {
		return nucl.pos + nucl.int_centers[nucl.BASE];
	}, R"pbdoc(
    Return the position of the base site.

    Returns
    -------
    numpy.ndarray
        The position of the site.
	)pbdoc");
}

#endif /* OXPY_BINDINGS_INCLUDES_BASEPARTICLE_H_ */
