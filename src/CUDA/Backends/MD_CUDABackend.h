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
#include "../Thermostats/CUDABrownianThermostat.h"
#include "../cuda_utils/cuda_device_utils.h"
#include "../Lists/CUDANoList.h"
#include "../Lists/CUDASimpleVerletList.h"

union CUDA_trap;

/**
 * @brief Manages a MD simulation on GPU with CUDA.
 */
class MD_CUDABackend: public MDBackend, public CUDABaseBackend{
protected:
	bool _use_edge;
	bool _any_rigid_body;
	bool _avoid_cpu_calculations;

	int *_h_gpu_index, *_h_cpu_index;

	std::shared_ptr<Timer> _timer_sorting;

	c_number4 *_d_vels, *_h_vels;
	c_number4 *_d_Ls, *_h_Ls;
	c_number4 *_d_forces, *_h_forces;
	c_number4 *_d_torques, *_h_torques;

	std::vector<int> _h_particles_to_mols;
	int *_d_particles_to_mols, *_d_mol_sizes, *_d_buff_particles_to_mols;
	c_number4 *_d_molecular_coms;

	c_number4 *_d_buff_vels, *_d_buff_Ls;

	llint _barostat_attempts, _barostat_accepted;
	int _update_st_every = 0;

	bool _print_energy;

	ObservableOutput *_obs_output_error_conf;
	std::string _error_conf_file;

	std::shared_ptr<CUDABaseThermostat> _cuda_thermostat;

	bool _cuda_barostat_always_refresh = false;
	std::shared_ptr<CUDABrownianThermostat> _cuda_barostat_thermostat;

	CUDA_trap *_d_ext_forces;
	int _max_ext_forces;

	virtual void _gpu_to_host();
	virtual void _host_to_gpu();
	virtual void _apply_external_forces_changes();

	virtual void _sort_particles();
	virtual void _rescale_molecular_positions(c_number4 new_Ls, c_number4 old_Ls, bool is_reverse_move);
	virtual void _rescale_positions(c_number4 new_Ls, c_number4 old_Ls);

	virtual void _first_step();
	virtual void _apply_barostat();
	virtual void _forces_second_step();
	virtual void _set_external_forces();

	virtual void _thermalize();
	virtual void _update_stress_tensor();

	virtual void _init_CUDA_MD_symbols();

public:
	MD_CUDABackend();
	virtual ~MD_CUDABackend();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void sim_step();

	virtual void apply_simulation_data_changes();
	virtual void apply_changes_to_simulation_data();
};

#endif /* MD_CUDABACKEND_H_ */
