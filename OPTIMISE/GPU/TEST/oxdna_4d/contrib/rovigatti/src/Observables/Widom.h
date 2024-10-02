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

class Widom: public BaseObservable {
protected:
	int _probe_type;
	Cells *_cells;
	BaseParticle *_probe;
	std::vector<BaseParticle *> _particles;
	int _N;
	int _tries;
	number _temperature;

public:
	Widom();
	virtual ~Widom();

	std::string get_output_string(llint curr_step);

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init();
};

extern "C" BaseObservable *make_Widom() {
	return new Widom();
}

#endif /* WIDOM_H_ */
