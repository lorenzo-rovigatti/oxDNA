/**
 * @file    MinBackend.h
 * @date    03/set/2010
 * @author  lorenzo
 *
 *
 */

#ifndef MINBACKEND_H_
#define MINBACKEND_H_

#include "MDBackend.h"
#include "./Thermostats/BaseThermostat.h"

/**
 * @brief Manages some kind of potential energy minimization
 */
template<typename number>
class MinBackend: public MDBackend<number> {
protected:
	void _compute_forces();
	void _evolve();

	number _max_step;

public:
	MinBackend();
	virtual ~MinBackend();

	void init();
	void get_settings (input_file &inp);
	void sim_step(llint cur_step);
	void activate_thermostat();
};

#endif /* MINBACKEND_H_ */
