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

class CUDANoThermostat: public CUDABaseThermostat, public NoThermostat {
public:
	CUDANoThermostat();
	virtual ~CUDANoThermostat();

	virtual void apply_cuda(c_number4 *d_poss, GPU_quat *d_orientationss, c_number4 *d_vels, c_number4 *d_Ls, llint curr_step);
	virtual bool would_activate(llint curr_step) {
		return false;
	}
};

#endif /* CUDANOTHERMOSTAT_H_ */
