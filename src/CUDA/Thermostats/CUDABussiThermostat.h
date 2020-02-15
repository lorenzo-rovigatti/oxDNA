/*
 * CUDABussiThermostat.h
 *
 *  Created on: MAy 09, 2015
 *      Author: rovigatti
 */

#ifndef CUDABUSSITHERMOSTAT_H_
#define CUDABUSSITHERMOSTAT_H_

#include "CUDABaseThermostat.h"
#include "../../Backends/Thermostats/BussiThermostat.h"
#include "../cuda_utils/CUDA_lr_common.cuh"

/**
 * @brief CUDA implementation of the {@link BussiThermostat thermostat} by Bussi et al.
 */

class CUDABussiThermostat: public CUDABaseThermostat, public BussiThermostat {
protected:
public:
	CUDABussiThermostat();
	virtual ~CUDABussiThermostat();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void apply_cuda(c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_vels, c_number4 *d_Ls, llint curr_step);
	virtual bool would_activate(llint curr_step);
};

#endif /* CUDABUSSITHERMOSTAT_H_ */
