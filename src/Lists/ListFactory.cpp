/*
 * ListFactory.cpp
 *
 *  Created on: 05/nov/2013
 *      Author: lorenzo
 */

#include "ListFactory.h"

#include "NoList.h"
#include "Cells.h"
#include "VerletList.h"
#include "BinVerletList.h"
#include "RodCells.h"

ListPtr ListFactory::make_list(input_file &inp, std::vector<BaseParticle *> &ps, BaseBox *box) {
	// the default list is verlet
	char list_type[512] = "verlet";
	getInputString(&inp, "list_type", list_type, 0);

	if(!strncmp(list_type, "verlet", 512)) return std::make_shared<VerletList>(ps, box);
	else if(!strncmp(list_type, "bin_verlet", 512)) return std::make_shared<BinVerletList>(ps, box);
	else if(!strncmp(list_type, "no", 512)) return std::make_shared<NoList>(ps, box);
	else if(!strncmp(list_type, "cells", 512)) return std::make_shared<Cells>(ps, box);
	else if(!strncmp(list_type, "rodcells", 512)) return std::make_shared<RodCells>(ps, box);
	else throw oxDNAException("Invalid list '%s'", list_type);
}
