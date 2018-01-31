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
#include "CUDATSPInteraction.h"
#include "CUDATEPInteraction.h"
#include "CUDARNAInteraction.h"

#include "../../Utilities/Utils.h"

CUDAInteractionFactory::CUDAInteractionFactory() {

}

CUDAInteractionFactory::~CUDAInteractionFactory() {

}

template<typename number, typename number4>
CUDABaseInteraction<number, number4> *CUDAInteractionFactory::make_interaction(input_file &inp) {
	// The default interaction is DNAInteraction
	string inter_type("DNA");
	getInputString(&inp, "interaction_type", inter_type, 0);

	if(!inter_type.compare("DNA") || !inter_type.compare("DNA_nomesh") || !inter_type.compare("DNA2")) return new CUDADNAInteraction<number, number4>();
	else if(!inter_type.compare("RNA") || !inter_type.compare("RNA2")  ) return new CUDARNAInteraction<number, number4>();
	else if(!inter_type.compare("LJ")) return new CUDALJInteraction<number, number4>();
	else if(!inter_type.compare("patchy")) return new CUDAPatchyInteraction<number, number4>();
	else if(!inter_type.compare("TSP")) return new CUDATSPInteraction<number, number4>();
	else if(inter_type.compare("TEP") == 0) return new CUDATEPInteraction<number, number4>();
	else {
		std::string cuda_name(inter_type);
		cuda_name = "CUDA" + cuda_name;
		CUDABaseInteraction<number, number4> *res = dynamic_cast<CUDABaseInteraction<number, number4> *>(PluginManager::instance()->get_interaction<number>(cuda_name));
		if(res == NULL) throw oxDNAException ("CUDA interaction '%s' not found. Aborting", cuda_name.c_str());
		return res;
	}
}

template class CUDABaseInteraction<float, float4> *CUDAInteractionFactory::make_interaction(input_file &inp);
template class CUDABaseInteraction<double, LR_double4> *CUDAInteractionFactory::make_interaction(input_file &inp);
