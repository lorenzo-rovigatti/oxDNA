/*
 * CUDAListFactory.cpp
 *
 *  Created on: 18/feb/2013
 *      Author: lorenzo
 */

#include "CUDAListFactory.h"
#include "../../Utilities/oxDNAException.h"

#include "CUDANoList.h"
#include "CUDASimpleVerletList.h"
#include "CUDABinVerletList.h"

CUDABaseList* CUDAListFactory::make_list(input_file &inp) {
	char list_type[256];

	if(getInputString(&inp, "CUDA_list", list_type, 0) == KEY_NOT_FOUND || !strcmp("verlet", list_type)) {
		return new CUDASimpleVerletList();
	}
	else if(!strcmp("no", list_type)) {
		return new CUDANoList();
	}
	else if(!strcmp("bin_verlet", list_type)) {
		return new CUDABinVerletList();
	}
	else {
		throw oxDNAException("CUDA_list '%s' is not supported", list_type);
	}
}
