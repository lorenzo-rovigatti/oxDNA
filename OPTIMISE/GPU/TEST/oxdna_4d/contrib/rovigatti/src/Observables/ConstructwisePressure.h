/*
 * ConstructwisePressure.h
 *
 *  Created on: 08/jan/2018
 *      Author: lorenzo
 */

#ifndef CONSTRUCTWISEPRESSURE_H_
#define CONSTRUCTWISEPRESSURE_H_

#include "Observables/BaseObservable.h"

/**
 * @brief Outputs the instantaneous pressure between constructs, computed through the virial.
 *
 * The supported syntax is
 * @verbatim
 type = pressure (an observable that computes the osmotic pressure of the system)
 construct_size = <int> (the size of the constructs)
 @endverbatim
 */

class ConstructwisePressure: public BaseObservable {
protected:
	double _P;
	number _shear_rate;
	int _construct_size;
	int _N_constructs;
	std::vector<LR_vector> _construct_coms;
	bool _is_polymer_swap = false;

public:
	ConstructwisePressure();
	virtual ~ConstructwisePressure();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	void init();

	void update_pressure();
	void update_pressure_PolymerSwap();
	
	number get_P() {
		return _P;
	}
	std::string get_output_string(llint curr_step);
};

extern "C" BaseObservable *make_ConstructwisePressure() {
	return new ConstructwisePressure();
}

#endif /* CONSTRUCTWISEPRESSURE_H_ */
