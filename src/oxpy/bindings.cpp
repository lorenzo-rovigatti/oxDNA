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
#include "bindings_includes/BaseBox.h"
#include "bindings_includes/Forces/BaseForce.h"
#include "bindings_includes/BaseInteraction.h"
#include "bindings_includes/Observables/BaseObservable.h"
#include "bindings_includes/BaseParticle.h"
#include "bindings_includes/FlattenedConfigInfo.h"
#include "bindings_includes/ConfigInfo.h"
#include "bindings_includes/DNANucleotide.h"
#include "bindings_includes/input_file.h"
#include "bindings_includes/RNANucleotide.h"
#include "bindings_includes/Molecule.h"

#include <Utilities/Utils.h>

PYBIND11_MODULE(core, m) {
	export_input_file(m);

	export_BaseBox(m);
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

	// export a utility function
	m.def("get_temperature", Utils::get_temperature, py::arg("raw_T"), R"pbdoc(
	Convert the given temperature string to its corresponding oxDNA value and return it.

	This function can be used to convert Celsius (*e.g.* `20 C`) or Kelvin (*e.g.* `300 K`) temperatures to oxDNA units. 

	Parameters
	----------
	raw_T: str
		The string storing the temperature.

	Returns
	-------
	float
		The temperature in oxDNA units.
	
)pbdoc");

	py::register_exception<oxDNAException>(m, "OxDNAError");

	py::module sub_m = m.def_submodule("analysis");
	export_AnalysisBackend(sub_m);
}
