/*
 * MetalInteractionFactory.mm
 *
 *  Created for Metal backend
 */

#include "MetalInteractionFactory.h"
#include "../../Utilities/oxDNAException.h"

#include "MetalDNAInteraction.h"
#include "MetalLJInteraction.h"
#include "MetalPatchyInteraction.h"
#include "MetalRNAInteraction.h"
#include "MetalTEPInteraction.h"

MetalBaseInteraction* MetalInteractionFactory::make_interaction(input_file &inp) {
	std::string interaction_type("DNA");
	getInputString(&inp, "interaction_type", interaction_type, 0);

	if(interaction_type == "DNA" || interaction_type == "DNA2" || interaction_type == "DNA_nomesh" || interaction_type == "DNA2_nomesh") {
		return new MetalDNAInteraction();
	}
	if(interaction_type == "LJ") return new MetalLJInteraction();
	if(interaction_type == "RNA" || interaction_type == "RNA2") return new MetalRNAInteraction();
	if(interaction_type == "patchy") return new MetalPatchyInteraction();
	if(interaction_type == "TEP") return new MetalTEPInteraction();

	throw oxDNAException("Metal interaction '%s' is not supported", interaction_type.c_str());
}
