/**
 * @file    FireBackend.h
 * @date    02/mar/2016
 * @author  flavio
 *
 *
 */

#ifndef FIREBACKEND_H_
#define FIREBACKEND_H_

#include "MDBackend.h"

/**
 * @brief Manages the FIRE algorithm for potential energy minimization
 */

class FIREBackend: public MDBackend {
protected:
	void _compute_forces();
	void _first_step();
	void _second_step();
	void _evolve();

	int _N_min, _step_last_reset;
	number _f_inc, _f_dec;
	number _alpha, _alpha_start, _f_alpha;
	number _max_step, _dt_max, _dt_restart;
	bool _allow_dt_increase;

public:
	FIREBackend();
	virtual ~FIREBackend();

	void init();
	void get_settings(input_file &inp);
	void sim_step();
	void activate_thermostat();
};

#endif /* FIREBACKEND_H_ */
