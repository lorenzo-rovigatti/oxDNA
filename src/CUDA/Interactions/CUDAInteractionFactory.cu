/*
 * CUDAInteractionFactory.cpp
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include <string>

#include "CUDAInteractionFactory.h"

#include "../../PluginManagement/PluginManager.h"

#include "CUDADNAInteraction.h"
#include "CUDALJInteraction.h"
#include "CUDAPatchyInteraction.h"
#include "CUDATEPInteraction.h"
#include "CUDARNAInteraction.h"

#include "../../Utilities/Utils.h"

CUDAInteractionFactory::CUDAInteractionFactory() {

}

CUDAInteractionFactory::~CUDAInteractionFactory() {

}

CUDABaseInteraction *CUDAInteractionFactory::make_interaction(input_file &inp) {
	// The default interaction is DNAInteraction
	string inter_type("DNA");
	getInputString(&inp, "interaction_type", inter_type, 0);

	if(!inter_type.compare("DNA") || !inter_type.compare("DNA_nomesh") || !inter_type.compare("DNA2")) return new CUDADNAInteraction();
	else if(!inter_type.compare("RNA") || !inter_type.compare("RNA2")  ) return new CUDARNAInteraction();
	else if(!inter_type.compare("LJ")) return new CUDALJInteraction();
	else if(!inter_type.compare("patchy")) return new CUDAPatchyInteraction();
	else if(inter_type.compare("TEP") == 0) return new CUDATEPInteraction();
	else {
		std::string cuda_name(inter_type);
		cuda_name = "CUDA" + cuda_name;
		CUDABaseInteraction *res = dynamic_cast<CUDABaseInteraction*>(PluginManager::instance()->get_interaction(cuda_name).get());
		if(res == NULL) throw oxDNAException ("CUDA interaction '%s' not found. Aborting", cuda_name.c_str());
		return res;
	}
}
