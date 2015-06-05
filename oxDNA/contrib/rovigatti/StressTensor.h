/*
 * StressTensor.h
 *
 *  Created on: 11/ott/2013
 *      Author: lorenzo
 */

#ifndef STRESSTENSOR_H_
#define STRESSTENSOR_H_

#include "BaseObservable.h"

/**
 * @brief Computes the stress tensors in real space (Eq. 4.4.13 in the Frenkel & Smit book).
 */
template<typename number>
class StressTensor: public BaseObservable<number> {
protected:
	int _N, _N_A, _N_B;
	bool _print_averaged_off_diagonal;

public:
	StressTensor();
	virtual ~StressTensor();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	std::string get_output_string(llint curr_step);
};

extern "C" BaseObservable<float> *make_float() { return new StressTensor<float>(); }
extern "C" BaseObservable<double> *make_double() { return new StressTensor<double>(); }

#endif /* STRESSTENSOR_H_ */
