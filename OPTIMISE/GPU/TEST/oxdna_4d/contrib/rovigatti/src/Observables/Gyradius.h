/*
 * Gyradius.h
 *
 *  Created on: Oct 30, 2013
 *      Author: Lorenzo
 */

#ifndef GYRADIUS_H_
#define GYRADIUS_H_

#include "Observables/BaseObservable.h"

/**
 * @brief Computes the gyration radius of a single chain
 *
 */

class Gyradius: public BaseObservable {
protected:
	bool _accumulate;
	number _avg_gr2;
	int _counter;
	int _type;

public:
	Gyradius();
	virtual ~Gyradius();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual std::string get_output_string(llint curr_step);
};

extern "C" BaseObservable *make_Gyradius() {
	return new Gyradius();
}

#endif /* GYRADIUS_H_ */
