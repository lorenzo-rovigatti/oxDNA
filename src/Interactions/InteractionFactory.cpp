/*
 * InteractionFactory.cpp
 *
 *  Created on: Feb 13, 2013
 *      Author: rovigatti
 */

#include "InteractionFactory.h"

#include "../PluginManagement/PluginManager.h"

#include "LJInteraction.h"
#include "DNAInteraction.h"
#include "DNA2Interaction.h"
#include "RNAInteraction.h"
#include "RNAInteraction2.h"
#include "DNAInteraction_nomesh.h"
#include "PatchyInteraction.h"
#include "PatchyInteractionDan.h"
#include "KFInteraction.h"
#include "DNAInteraction_relax.h"
#include "HSInteraction.h"
#include "DHSInteraction.h"
#include "BoxInteraction.h"
#include "HardCylinderInteraction.h"
#include "HardSpheroCylinderInteraction.h"
#include "CustomInteraction.h"
#include "RNAInteraction_relax.h"
#include "TEPInteraction.h"
#include "JordanInteraction.h"
#include "DRHInteraction.h"
#include "DRHInteraction_relax.h"


InteractionPtr InteractionFactory::make_interaction(input_file &inp) {
	// The default interaction is DNAInteraction
	std::string inter_type("DNA");
	getInputString(&inp, "interaction_type", inter_type, 0);
	std::string backend;
	getInputString(&inp, "backend", backend, 0);

	if(inter_type.compare("DNA") == 0) {
		// in order to avoid small mismatches between potential energies computed on the GPU and 
		// on the CPU we enforce the DNAInteraction_nomesh interaction class
		if(backend.compare("CUDA") == 0) return std::make_shared<DNAInteraction_nomesh>();
		else return std::make_shared<DNAInteraction>();
	}
	else if(inter_type.compare("DNA2") == 0) {
		if(backend.compare("CUDA") == 0) return std::make_shared<DNA2Interaction_nomesh>();
		else return std::make_shared<DNA2Interaction>();
	}
	else if(inter_type.compare("LJ") == 0) return std::make_shared<LJInteraction>();
	else if(inter_type.compare("DNA_nomesh") == 0) return std::make_shared<DNAInteraction_nomesh>();
	else if(inter_type.compare("DNA2_nomesh") == 0) return std::make_shared<DNA2Interaction_nomesh>();
	else if(inter_type.compare("DNA_relax") == 0) return std::make_shared<DNAInteraction_relax>();
	else if(inter_type.compare("DNA2") == 0) return std::make_shared<DNA2Interaction>();
	else if(inter_type.compare("RNA") == 0) return std::make_shared<RNAInteraction>();
	else if(inter_type.compare("RNA2") == 0) return std::make_shared<RNA2Interaction>();
	else if(inter_type.compare("RNA_relax") == 0) return std::make_shared<RNAInteraction_relax>();
	else if(inter_type.compare("patchy") == 0) return std::make_shared<PatchyInteraction>();
	else if(inter_type.compare("patchyDan") == 0) return std::make_shared<PatchyInteractionDan>();
	else if(inter_type.compare("KF") == 0) return std::make_shared<KFInteraction>();
	else if(inter_type.compare("HS") == 0) return std::make_shared<HSInteraction>();
	else if(inter_type.compare("Box") == 0) return std::make_shared<BoxInteraction>();
	else if(inter_type.compare("HardCylinder") == 0) return std::make_shared<HardCylinderInteraction>();
	else if(inter_type.compare("HardSpheroCylinder") == 0) return std::make_shared<HardSpheroCylinderInteraction>();
	else if(inter_type.compare("DHS") == 0) return std::make_shared<DHSInteraction>();
	else if(inter_type.compare("custom") == 0) return std::make_shared<CustomInteraction>();
	else if(inter_type.compare("TEP") == 0) return std::make_shared<TEPInteraction>();
	else if(inter_type.compare("Jordan") == 0) return std::make_shared<JordanInteraction>();
	else if(inter_type.compare("NA") == 0) return std::make_shared<DRHInteraction>();
	else if(inter_type.compare("NA_relax") == 0) return std::make_shared<DRHInteraction_relax>();
	else {
		InteractionPtr res = PluginManager::instance()->get_interaction(inter_type);
		if(res == NULL) throw oxDNAException("Interaction '%s' not found. Aborting", inter_type.c_str());
		return res;
	}
}
