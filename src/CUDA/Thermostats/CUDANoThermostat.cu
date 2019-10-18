/*
 * CUDANoThermostat.cu
 *
 *  Created on: Feb 15, 2013
 *      Author: rovigatti
 */

#include "CUDANoThermostat.h"

CUDANoThermostat::CUDANoThermostat() :
				NoThermostat() {

}

CUDANoThermostat::~CUDANoThermostat() {

}

void CUDANoThermostat::apply_cuda(number4 *d_poss, GPU_quat *d_orientationss, number4 *d_vels, number4 *d_Ls, llint curr_step) {
	return;
}
