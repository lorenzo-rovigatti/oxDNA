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
#include "./Thermostats/BaseThermostat.h"

/**
 * @brief Manages a MD simulation on CPU. It supports NVE and NVT (Brownian or Langevin) simulations
 */
template<typename number>
class MD_CPUBackend: public MDBackend<number> {
protected:
	LR_vector<number> _vcm;

	BaseThermostat<number> * _thermostat;

	void _first_step(llint cur_step);
	void _compute_forces();
	void _second_step();

public:
	MD_CPUBackend();
	virtual ~MD_CPUBackend();

	void init();
	void get_settings (input_file &inp);
	void sim_step(llint cur_step);
	void activate_thermostat();
};

#endif /* MD_CPUBACKEND_H_ */
