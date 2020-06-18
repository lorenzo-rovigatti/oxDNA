/*
 * CUDABrownianThermostat.h
 *
 *  Created on: Feb 15, 2013
 *      Author: rovigatti
 */

#ifndef CUDABROWNIANTHERMOSTAT_H_
#define CUDABROWNIANTHERMOSTAT_H_

#include "CUDABaseThermostat.h"
#include "../../Backends/Thermostats/BrownianThermostat.h"

/**
 * @brief CUDA implementation of a {@link BrownianThermostat brownian thermostat}.
 */

class CUDABrownianThermostat: public CUDABaseThermostat, public BrownianThermostat {
protected:
public:
	CUDABrownianThermostat();
	virtual ~CUDABrownianThermostat();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void apply_cuda(c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_vels, c_number4 *d_Ls, llint curr_step);
	virtual bool would_activate(llint curr_step);
};

#endif /* CUDABROWNIANTHERMOSTAT_H_ */
