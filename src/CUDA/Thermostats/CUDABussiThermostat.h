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
template<typename number, typename number4>
class CUDABussiThermostat: public CUDABaseThermostat<number, number4>, public BussiThermostat<number> {
protected:
public:
	CUDABussiThermostat();
	virtual ~CUDABussiThermostat();

	virtual void get_settings(input_file &inp);
	virtual void init(int N);

	virtual void apply_cuda(number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_vels, number4 *d_Ls, llint curr_step);
	virtual bool would_activate(llint curr_step);
};

#endif /* CUDABUSSITHERMOSTAT_H_ */
