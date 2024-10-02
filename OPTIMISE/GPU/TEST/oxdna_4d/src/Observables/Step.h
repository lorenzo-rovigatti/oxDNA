/*
 * Step.h
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#ifndef STEP_H_
#define STEP_H_

#include "BaseObservable.h"

/**
 * @brief Outputs the simulation step (optionally scaling it by dt).
 *
 * This observable takes one optional argument:
 *
 * @verbatim
 [units = steps|MD (units to print the time on. time in MD units = steps * dt, defaults to step)]
 @endverbatim
 */

class Step: public BaseObservable {
private:
	enum {
		STEP_UNIT_ONE_PER_STEP = 1, STEP_UNIT_HONOUR_DT = 2
	};

	number _dt;
	int _units;
public:
	Step();
	virtual ~Step();

	virtual std::string get_output_string(llint curr_step);
	void get_settings(input_file &my_inp, input_file &sim_inp);
};

#endif /* STEP_H_ */
