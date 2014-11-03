/*
 * CUDALangevinThermostat.h
 *
 *  Created on: Sep 1, 2014
 *      Author: mzimmer
 */

#ifndef CUDALANGEVINTHERMOSTAT_H_
#define CUDALANGEVINTHERMOSTAT_H_

#include "CUDABaseThermostat.h"
#include "../../Backends/Thermostats/LangevinThermostat.h"
#include "../cuda_utils/cuda_device_utils.h"
#include "../cuda_utils/CUDA_lr_common.cuh"

template<typename number, typename number4>
class CUDALangevinThermostat: public CUDABaseThermostat<number, number4>, public LangevinThermostat<number> {
protected:
public:
	CUDALangevinThermostat();
	virtual ~CUDALangevinThermostat();

	virtual void get_settings(input_file &inp);
	virtual void init(int N);

	virtual void apply_cuda(number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_vels, number4 *d_Ls, llint curr_step);
	virtual bool would_activate(llint curr_step);
};

#endif /* CUDALANGEVINTHERMOSTAT_H_ */
