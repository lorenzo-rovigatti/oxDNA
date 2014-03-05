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
[stress_tensor = <bool> (if true, the output will contain 7 fields, with the first being the total pressure and the other 6 the six independent components of the stress tensor, xx, yy, zz, xy, xz, yz)]
@endverbatim
 */
template<typename number>
class Pressure: public BaseObservable<number> {
protected:
	number _T;
	bool _stress_tensor;

public:
	Pressure();
	virtual ~Pressure();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	std::string get_output_string(llint curr_step);
};

#endif /* PRESSURE_H_ */
