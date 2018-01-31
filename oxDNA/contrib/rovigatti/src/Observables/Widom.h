/*
 * Widom.h
 *
 *  Created on: Jan 14, 2016
 *      Author: Lorenzo Rovigatti
 */

#ifndef WIDOM_H_
#define WIDOM_H_

#include "Lists/Cells.h"

#include "Observables/BaseObservable.h"
#include "Interactions/DNAInteraction.h"
#include "Utilities/OrderParameters.h"

template<typename number>
class Widom : public BaseObservable<number>  {
protected:
	int _probe_type;
	Cells<number> *_cells;
	BaseParticle<number> *_probe;
	BaseParticle<number> **_particles;
	int _N;
	int _tries;
	number _temperature;

public:
	Widom();
	virtual ~Widom();

	std::string get_output_string(llint curr_step);

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init(ConfigInfo<number> &config_info);
};

extern "C" BaseObservable<float> *make_float() { return new Widom<float>(); }
extern "C" BaseObservable<double> *make_double() { return new Widom<double>(); }

#endif /* WIDOM_H_ */
