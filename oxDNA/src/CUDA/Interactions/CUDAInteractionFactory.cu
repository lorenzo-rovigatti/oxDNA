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
#include "../../Utilities/Utils.h"
#include "CUDARNAInteraction.h"


CUDAInteractionFactory::CUDAInteractionFactory() {

}

CUDAInteractionFactory::~CUDAInteractionFactory() {

}

template<typename number, typename number4>
CUDABaseInteraction<number, number4> *CUDAInteractionFactory::make_interaction(input_file &inp) {
	// The default interaction is DNAInteraction
	char inter_type[512] = "DNA";
	getInputString(&inp, "interaction_type", inter_type, 0);

	if(!strncmp(inter_type, "DNA", 512) || !strncmp(inter_type, "DNA_nomesh", 512) || !strncmp(inter_type, "DNA2", 512)) return new CUDADNAInteraction<number, number4>();
	else if(!strncmp(inter_type, "RNA", 512) || !strncmp(inter_type, "RNA2", 512)  ) return new CUDARNAInteraction<number, number4>();
	else if(!strncmp(inter_type, "LJ", 512)) return new CUDALJInteraction<number, number4>();
	else if(!strncmp(inter_type, "patchy", 512)) return new CUDAPatchyInteraction<number, number4>();
	else if(!strncmp(inter_type, "TSP", 512)) return new CUDATSPInteraction<number, number4>();
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
