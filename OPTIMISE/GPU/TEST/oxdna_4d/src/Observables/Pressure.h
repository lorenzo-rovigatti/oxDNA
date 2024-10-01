/*
 * Pressure.h
 *
 *  Created on: 25/ott/2013
 *      Author: lorenzo
 */

#ifndef PRESSURE_H_
#define PRESSURE_H_

#include "BaseObservable.h"

/**
 * @brief Outputs the instantaneous pressure, computed through the virial.
 *
 * The supported syntax is
 * @verbatim
 type = pressure (an observable that computes the osmotic pressure of the system)
 [stress_tensor = <bool> (if true, the output will contain 9 fields: the total pressure and the nine components of the stress tensor, xx, xy, xz, yx, yy, yz, zx, zy, zz)]
 @endverbatim
 */

class Pressure: public BaseObservable {
protected:
	bool _custom_stress_tensor;
	bool _with_stress_tensor;
	bool _PV_only;

public:
	Pressure();
	virtual ~Pressure();

	void get_settings(input_file &my_inp, input_file &sim_inp);

	std::string get_output_string(llint curr_step);
};

#endif /* PRESSURE_H_ */
