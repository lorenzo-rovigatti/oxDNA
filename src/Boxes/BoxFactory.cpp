/*
 * BoxFactory.cpp
 *
 *  Created on: 17/mar/2013
 *      Author: lorenzo
 */

#include "BoxFactory.h"

#include "CubicBox.h"
#include "OrthogonalBox.h"
#include "LeesEdwardsCubicBox.h"

#include "../Utilities/oxDNAException.h"

BoxPtr BoxFactory::make_box(input_file &inp) {
	// the default box is the cubic one
	char box_type[512] = "cubic";
	getInputString(&inp, "box_type", box_type, 0);
	bool lees_edwards = false;
	getInputBool(&inp, "lees_edwards", &lees_edwards, 0);

	if(!strncmp(box_type, "cubic", 512)) {
		if(!lees_edwards) return std::make_shared<CubicBox>();
		else return std::make_shared<LeesEdwardsCubicBox>();
	}
	else if(!strncmp(box_type, "orthogonal", 512)) {
		if(!lees_edwards) return std::make_shared<OrthogonalBox>();
		else throw oxDNAException("Lees-Edwards boundary conditions and orthogonal boxes are not compatible");
	}
	else throw oxDNAException("Unsupported box '%s'", box_type);
}
