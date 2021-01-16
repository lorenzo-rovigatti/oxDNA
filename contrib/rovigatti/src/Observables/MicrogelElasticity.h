/*
 * MicrogelElasticity.h
 *
 *  Created on: Mar 13, 2018
 *      Author: Lorenzo
 */

#ifndef MICROGELELASTICITY_H_
#define MICROGELELASTICITY_H_

#include "Observables/BaseObservable.h"

class MicrogelElasticity: public BaseObservable {
protected:
	LR_vector _com();
	bool _crosslinkers_only;
	int _N_crosslinkers;

public:
	MicrogelElasticity();
	virtual ~MicrogelElasticity();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	void init();
	std::string get_output_string(llint curr_step);
};

extern "C" BaseObservable *make_MicrogelElasticity() {
	return new MicrogelElasticity();
}

#endif /* MICROGELELASTICITY_H_ */
