/*
 * CUDANoThermostat.h
 *
 *  Created on: Feb 15, 2013
 *      Author: rovigatti
 */

#ifndef CUDANOTHERMOSTAT_H_
#define CUDANOTHERMOSTAT_H_

#include "CUDABaseThermostat.h"
#include "../../Backends/Thermostats/NoThermostat.h"

/**
 * @brief Implements a NVE simulation (no thermalisation).
 */
template<typename number, typename number4>
class CUDANoThermostat: public CUDABaseThermostat<number, number4>, public NoThermostat<number> {
public:
	CUDANoThermostat();
	virtual ~CUDANoThermostat();

	virtual void apply_cuda(number4 *d_poss, GPU_quat<number> *d_orientationss, number4 *d_vels, number4 *d_Ls, llint curr_step);
	virtual bool would_activate(llint curr_step) { return false; }
};

#endif /* CUDANOTHERMOSTAT_H_ */
