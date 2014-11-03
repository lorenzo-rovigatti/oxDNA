/*
 * CUDANoThermostat.cu
 *
 *  Created on: Feb 15, 2013
 *      Author: rovigatti
 */

#include "CUDANoThermostat.h"

template<typename number, typename number4>
CUDANoThermostat<number, number4>::CUDANoThermostat() : NoThermostat<number>() {

}

template<typename number, typename number4>
CUDANoThermostat<number, number4>::~CUDANoThermostat() {

}

template<typename number, typename number4>
void CUDANoThermostat<number, number4>::apply_cuda(number4 *d_poss, GPU_quat<number> *d_orientationss, number4 *d_vels, number4 *d_Ls, llint curr_step) {
	return;
}

template class CUDANoThermostat<float, float4>;
template class CUDANoThermostat<double, LR_double4>;
