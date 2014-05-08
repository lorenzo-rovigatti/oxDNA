/*
 * Gyradius.h
 *
 *  Created on: Oct 30, 2013
 *      Author: Lorenzo
 */

#ifndef GYRADIUS_H_
#define GYRADIUS_H_

#include "BaseObservable.h"

/**
 * @brief Computes the gyration radius of a single chain
 *
 */
template<typename number>
class Gyradius : public BaseObservable<number> {
protected:
	bool _accumulate;
	number _avg_gr2;
	int _counter;
	int _type;

public:
	Gyradius();
	virtual ~Gyradius();

	void get_settings (input_file &my_inp, input_file &sim_inp);
	virtual std::string get_output_string(llint curr_step);
};

extern "C" BaseObservable<float> *make_float() { return new Gyradius<float>(); }
extern "C" BaseObservable<double> *make_double() { return new Gyradius<double>(); }

#endif /* GYRADIUS_H_ */
