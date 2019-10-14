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

/*#include "DNAInteraction_nomesh.h"
#include "RNAInteraction.h"
#include "PatchyInteraction.h"
#include "PatchyInteractionDan.h"
#include "KFInteraction.h"
#include "DNAInteraction_relax.h"
#include "TSPInteraction.h"
#include "HSInteraction.h"
#include "DHSInteraction.h"
#include "BoxInteraction.h"
#include "HardCylinderInteraction.h"
#include "HardSpheroCylinderInteraction.h"
#include "CustomInteraction.h"
#include "RNAInteraction2.h"
#include "RNAInteraction_relax.h"
#include "TEPInteraction.h"
#include "JordanInteraction.h"
*/

InteractionFactory::InteractionFactory() {

}

InteractionFactory::~InteractionFactory() {

}


IBaseInteraction *InteractionFactory::make_interaction(input_file &inp) {
	// The default interaction is DNAInteraction
	std::string inter_type ("DNA");
	getInputString(&inp, "interaction_type", inter_type, 0);

	if(inter_type.compare("DNA") == 0) {
		// in order to avoid small mismatches between potential energies computed on the GPU and 
		// on the CPU we enforce the DNAInteraction_nomesh interaction class
		std::string backend ("");
		getInputString(&inp, "backend", backend, 0);
		/*if(backend.compare("CUDA") == 0) return new DNAInteraction_nomesh();
		else */return new DNAInteraction();
	}
	else if(inter_type.compare("DNA2") == 0) {
		std::string backend ("");
		getInputString(&inp, "backend", backend, 0);
		/*if(backend.compare("CUDA") == 0) return new DNA2Interaction_nomesh();
		else */return new DNA2Interaction();
	}
	else if(inter_type.compare("LJ") == 0) return new LJInteraction();
	/*else if(inter_type.compare("DNA_nomesh") == 0) return new DNAInteraction_nomesh();
	else if(inter_type.compare("DNA2_nomesh") == 0) return new DNA2Interaction_nomesh();
	else if(inter_type.compare("DNA_relax") == 0) return new DNAInteraction_relax();
	else if(inter_type.compare("RNA") == 0) return new RNAInteraction();
	else if(inter_type.compare("patchy") == 0) return new PatchyInteraction();
	else if(inter_type.compare("patchyDan") == 0) return new PatchyInteractionDan();
	else if(inter_type.compare("KF") == 0) return new KFInteraction();
	else if(inter_type.compare("TSP") == 0) return new TSPInteraction();
	else if(inter_type.compare("HS") == 0) return new HSInteraction();
	else if(inter_type.compare("Box") == 0) return new BoxInteraction();
	else if(inter_type.compare("HardCylinder") == 0) return new HardCylinderInteraction();
	else if(inter_type.compare("HardSpheroCylinder") == 0) return new HardSpheroCylinderInteraction();
	else if(inter_type.compare("DHS") == 0) return new DHSInteraction();
	else if(inter_type.compare("custom") == 0) return new CustomInteraction();
	else if(inter_type.compare("DNA2") == 0) return new DNA2Interaction();
	else if(inter_type.compare("RNA2") == 0) return new RNA2Interaction();
	else if(inter_type.compare("RNA_relax") == 0) return new RNAInteraction_relax();
	else if(inter_type.compare("TEP") == 0) return new TEPInteraction();
	else if(inter_type.compare("Jordan") == 0) return new JordanInteraction();*/
	else {
		IBaseInteraction *res = PluginManager::instance()->get_interaction(inter_type);
		if(res == NULL) throw oxDNAException ("Interaction '%s' not found. Aborting", inter_type.c_str());
		return res;
	}
}
