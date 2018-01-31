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
template<typename number>
class ConstructwisePressure: public BaseObservable<number> {
protected:
	number _T;
	double _P;
	number _shear_rate;
	int _construct_size;
	int _N_constructs;
	vector<LR_vector<number> > _construct_coms;

public:
	ConstructwisePressure();
	virtual ~ConstructwisePressure();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	void init(ConfigInfo<number> &Info);

	void update_pressure();
	number get_P() { return _P; }
	std::string get_output_string(llint curr_step);
};

extern "C" BaseObservable<float> *make_float() { return new ConstructwisePressure<float>(); }
extern "C" BaseObservable<double> *make_double() { return new ConstructwisePressure<double>(); }

#endif /* CONSTRUCTWISEPRESSURE_H_ */
