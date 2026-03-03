/*
 * CUDARefreshThermostat.h
 *
 *  CUDA implementation of RefreshThermostat
 */

#ifndef CUDAREFRESHTHERMOSTAT_H_
#define CUDAREFRESHTHERMOSTAT_H_

#include "CUDABaseThermostat.h"
#include "../../Backends/Thermostats/RefreshThermostat.h"

/**
 * @brief CUDA implementation of a {@link RefreshThermostat refresh thermostat}.
 *
 * Refreshes every particle's translational and rotational momenta every n timesteps.
 */
class CUDARefreshThermostat: public CUDABaseThermostat, public RefreshThermostat {
public:
	CUDARefreshThermostat();
	virtual ~CUDARefreshThermostat();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void apply_cuda(c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_vels, c_number4 *d_Ls, llint curr_step);
	virtual bool would_activate(llint curr_step);
};

#endif /* CUDAREFRESHTHERMOSTAT_H_ */