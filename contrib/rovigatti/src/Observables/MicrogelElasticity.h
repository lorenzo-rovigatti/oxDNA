/*
 * MicrogelElasticity.h
 *
 *  Created on: Mar 13, 2018
 *      Author: Lorenzo
 */

#ifndef MICROGELELASTICITY_H_
#define MICROGELELASTICITY_H_

#include "Observables/BaseObservable.h"

template<typename number>
class MicrogelElasticity : public BaseObservable<number> {
protected:
	LR_vector<number> _com();
	bool _crosslinkers_only;
	int _N_crosslinkers;

public:
	MicrogelElasticity();
	virtual ~MicrogelElasticity();

	void get_settings (input_file &my_inp, input_file &sim_inp);
	void init(ConfigInfo<number> &config_info);
	std::string get_output_string(llint curr_step);
};

extern "C" BaseObservable<float> *make_MicrogelElasticity_float() { return new MicrogelElasticity<float>(); }
extern "C" BaseObservable<double> *make_MicrogelElasticity_double() { return new MicrogelElasticity<double>(); }

#endif /* MICROGELELASTICITY_H_ */
