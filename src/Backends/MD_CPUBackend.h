/**
 * @file    MD_CPUBackend.h
 * @date    03/set/2010
 * @author  lorenzo
 *
 *
 */

#ifndef MD_CPUBACKEND_H_
#define MD_CPUBACKEND_H_

#include "MDBackend.h"
#include "MCMoves/VolumeMove.h"

class BaseThermostat;

/**
 * @brief Manages a MD simulation on CPU. It supports NVE and NVT simulations
 */

class MD_CPUBackend: public MDBackend {
protected:
	std::shared_ptr<BaseThermostat> _thermostat;
	MovePtr _V_move;
	int _stress_tensor_avg_every;
	int _stress_tensor_counter;

	// thermostat introduced in https://journals.aps.org/pre/abstract/10.1103/PhysRevE.75.056707
	bool _use_builtin_langevin_thermostat = false;
	number _langevin_c1 = 0.;
	number _langevin_c2 = 0.;

	void _first_step();
	void _compute_forces();
	void _second_step();

	void _update_backend_info();

public:
	MD_CPUBackend();
	virtual ~MD_CPUBackend();

	void init();
	void get_settings(input_file &inp);
	void sim_step();
	void activate_thermostat();
};

#endif /* MD_CPUBACKEND_H_ */
