/*
 * MetalNoList.mm
 *
 *  Created for Metal backend
 */

#include "MetalNoList.h"

#include "../../Utilities/oxDNAException.h"

MetalNoList::MetalNoList() {

}

MetalNoList::~MetalNoList() {

}

void MetalNoList::get_settings(input_file &inp) {
	bool use_edge = false;
	getInputBool(&inp, "use_edge", &use_edge, 0);
	if(use_edge) {
		throw oxDNAException("'Metal_list = no' and 'use_edge = true' are incompatible");
	}
}
