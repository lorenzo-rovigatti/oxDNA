/*
 * CUDANoList.cpp
 *
 *  Created on: 29/set/2010
 *      Author: lorenzo
 */

#include "CUDANoList.h"

#include "../../Utilities/oxDNAException.h"

CUDANoList::CUDANoList() {

}

CUDANoList::~CUDANoList() {

}

void CUDANoList::get_settings(input_file &inp) {
	bool use_edge = false;
	getInputBool(&inp, "use_edge", &use_edge, 0);
	if(use_edge) {
		throw oxDNAException("'CUDA_list = no' and 'use_edge = true' are incompatible");
	}
}
