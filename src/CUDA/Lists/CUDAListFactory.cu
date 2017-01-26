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

CUDAListFactory::CUDAListFactory() {

}

CUDAListFactory::~CUDAListFactory() {

}

template <typename number, typename number4>
CUDABaseList<number, number4> *CUDAListFactory::make_list(input_file &inp) {
	char list_type[256];

	if(getInputString(&inp, "CUDA_list", list_type, 0) == KEY_NOT_FOUND || !strcmp("no", list_type)) return new CUDANoList<number, number4>();
	else if(!strcmp("verlet", list_type)) return new CUDASimpleVerletList<number, number4>();
	else if(!strcmp("bin_verlet", list_type)) return new CUDABinVerletList<number, number4>();
	else throw new oxDNAException("CUDA_list '%s' is not supported", list_type);
}

template CUDABaseList<float, float4> *CUDAListFactory::make_list(input_file &inp);
template CUDABaseList<double, LR_double4> *CUDAListFactory::make_list(input_file &inp);
