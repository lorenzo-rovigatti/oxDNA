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

void CUDANoThermostat::apply_cuda(tmpnmbr *d_poss, GPU_quat *d_orientationss, tmpnmbr *d_vels, tmpnmbr *d_Ls, llint curr_step) {
	return;
}
