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
#include "../cuda_utils/CUDA_lr_common.cuh"

/**
 * @brief CUDA implementation of a {@link BrownianThermostat brownian thermostat}.
 */
template<typename number, typename number4>
class CUDABrownianThermostat: public CUDABaseThermostat<number, number4>, public BrownianThermostat<number> {
protected:
public:
	CUDABrownianThermostat();
	virtual ~CUDABrownianThermostat();

	virtual void get_settings(input_file &inp);
	virtual void init(int N);

	virtual void apply_cuda(number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_vels, number4 *d_Ls, llint curr_step);
	virtual bool would_activate(llint curr_step);
};

#endif /* CUDABROWNIANTHERMOSTAT_H_ */
