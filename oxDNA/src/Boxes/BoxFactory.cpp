/*
 * BoxFactory.cpp
 *
 *  Created on: 17/mar/2013
 *      Author: lorenzo
 */

#include "BoxFactory.h"
#include "CubicBox.h"
#include "OrthogonalBox.h"

BoxFactory::BoxFactory() {

}

BoxFactory::~BoxFactory() {

}

template<typename number>
BaseBox<number> *BoxFactory::make_box(input_file &inp) {
	// the default box is the cubic one
	char box_type[512] = "cubic";
	getInputString(&inp, "box_type", box_type, 0);

	if(!strncmp(box_type, "cubic", 512)) return new CubicBox<number>();
	else if(!strncmp(box_type, "orthogonal", 512)) return new OrthogonalBox<number>();
	else throw oxDNAException("Unsupported box '%s'", box_type);
}

template BaseBox<float> *BoxFactory::make_box<float>(input_file &inp);
template BaseBox<double> *BoxFactory::make_box<double>(input_file &inp);
