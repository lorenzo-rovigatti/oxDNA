/*
 * bindings.cpp
 *
 *  Created on: Sep 13, 2019
 *      Author: lorenzo
 */

#include "python_defs.h"

#include "OxpyContext.h"
#include "OxpyManager.h"

#include "vector_matrix_casters.h"

#include "bindings_includes/AnalysisBackend.h"
#include "bindings_includes/BaseForce.h"
#include "bindings_includes/BaseInteraction.h"
#include "bindings_includes/BaseObservable.h"
#include "bindings_includes/BaseParticle.h"
#include "bindings_includes/FlattenedConfigInfo.h"
#include "bindings_includes/ConfigInfo.h"
#include "bindings_includes/DNANucleotide.h"
#include "bindings_includes/input_file.h"
#include "bindings_includes/RNANucleotide.h"
#include "bindings_includes/Molecule.h"

PYBIND11_MODULE(core, m) {
	export_input_file(m);

	export_BaseObservable(m);
	export_BaseParticle(m);
	export_BaseForce(m);
	export_FlattenedConfigInfo(m);
	export_ConfigInfo(m);
	export_DNANucleotide(m);
	export_BaseInteraction(m);
	export_RNANucleotide(m);

	export_Molecule(m);

	export_OxpyContext(m);

	export_SimManager(m);
	export_OxpyManager(m);

	py::register_exception<oxDNAException>(m, "OxDNAError");

	py::module sub_m = m.def_submodule("analysis");
	export_AnalysisBackend(sub_m);
}
