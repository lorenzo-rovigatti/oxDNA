/**
 * @file    MD_CUDABackend.h
 * @date    03/set/2010
 * @author  lorenzo
 *
 *
 */

#ifndef MD_CUDABACKEND_H_
#define MD_CUDABACKEND_H_

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "CUDABaseBackend.h"
#include "../../Backends/MDBackend.h"

#include "../CUDAUtils.h"
#include "../Thermostats/CUDABaseThermostat.h"
#include "../cuda_utils/cuda_device_utils.h"
#include "../Lists/CUDANoList.h"
#include "../Lists/CUDASimpleVerletList.h"

#include "../CUDAForces.h"

/**
 * @brief Manages a MD simulation on GPU with CUDA.
 */

class MD_CUDABackend: public MDBackend, public CUDABaseBackend{
protected:
	bool _use_edge;
	bool _any_rigid_body;
	bool _restart_step_counter;
	bool _avoid_cpu_calculations;

	int *_h_gpu_index, *_h_cpu_index;

	std::shared_ptr<Timer> _timer_sorting;

	c_number4 *_d_vels, *_h_vels;
	c_number4 *_d_Ls, *_h_Ls;
	c_number4 *_d_forces, *_h_forces;
	c_number4 *_d_torques, *_h_torques;

	c_number4 *_d_buff_vels, *_d_buff_Ls;
	llint _curr_step;

	llint _barostat_attempts, _barostat_accepted;

	bool _print_energy;

	ObservableOutput *_obs_output_error_conf;
	std::string _error_conf_file;

	CUDABaseThermostat*_cuda_thermostat;

	//constant_rate_force *_h_ext_forces, *_d_ext_forces;
	//mutual_trap *_h_ext_forces, *_d_ext_forces;
	CUDA_trap *_h_ext_forces, *_d_ext_forces;
	int _max_ext_forces;

	virtual void _host_to_gpu();
	virtual void _gpu_to_host();

	virtual void _host_particles_to_gpu();
	virtual void _gpu_to_host_particles();

	virtual void _sort_particles();
	virtual void _rescale_positions(c_number4 new_Ls, c_number4 old_Ls);

	virtual void _first_step();
	virtual void _apply_barostat(llint curr_step);
	virtual void _forces_second_step();
	virtual void _set_external_forces();

	virtual void _thermalize(llint curr_step);

	virtual void _init_CUDA_MD_symbols();

	virtual void _print_ready_observables(llint curr_step);

public:
	MD_CUDABackend();
	virtual ~MD_CUDABackend();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void sim_step(llint curr_step);
	virtual void fix_diffusion();

	virtual void print_conf(llint curr_step, bool reduced=false, bool only_last=false);
};

#endif /* MD_CUDABACKEND_H_ */
