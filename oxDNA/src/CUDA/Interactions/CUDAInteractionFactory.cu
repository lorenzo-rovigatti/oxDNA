/*
 * CUDAInteractionFactory.cpp
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include "CUDAInteractionFactory.h"

#include "CUDADNAInteraction.h"
#include "CUDALJInteraction.h"
#include "CUDAPatchyInteraction.h"
#include "CUDATSPInteraction.h"
#include "../../Utilities/Utils.h"

CUDAInteractionFactory::CUDAInteractionFactory() {

}

CUDAInteractionFactory::~CUDAInteractionFactory() {

}

template<typename number, typename number4>
CUDABaseInteraction<number, number4> *CUDAInteractionFactory::make_interaction(input_file &inp) {
	// The default interaction is DNAInteraction
	char inter_type[512] = "DNA";
	getInputString(&inp, "interaction_type", inter_type, 0);

	if(!strncmp(inter_type, "DNA", 512) || !strncmp(inter_type, "DNA_nomesh", 512)) return new CUDADNAInteraction<number, number4>();
	else if(!strncmp(inter_type, "LJ", 512)) return new CUDALJInteraction<number, number4>();
	else if(!strncmp(inter_type, "patchy", 512)) return new CUDAPatchyInteraction<number, number4>();
	else if(!strncmp(inter_type, "TSP", 512)) return new CUDATSPInteraction<number, number4>();
	else throw oxDNAException("Invalid interaction '%s'", inter_type);
}

template class CUDABaseInteraction<float, float4> *CUDAInteractionFactory::make_interaction(input_file &inp);
template class CUDABaseInteraction<double, LR_double4> *CUDAInteractionFactory::make_interaction(input_file &inp);
