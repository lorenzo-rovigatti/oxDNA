/*
 * Remoteness.h
 *
 *  Created on: Sep 23, 2016
 *      Author: Lorenzo Rovigatti
 */

#ifndef REMOTENESS_H_
#define REMOTENESS_H_

#include "Lists/Cells.h"

#include "Observables/BaseObservable.h"
#include "Interactions/DNAInteraction.h"
#include "Utilities/OrderParameters.h"

class Remoteness: public BaseObservable {
protected:
	Cells *_cells;
	BaseParticle *_probe;
	std::vector<BaseParticle *> _particles;
	int _N;
	int _tries;
	int _n_bins;
	int _n_confs;
	number _bin_size;
	number _max_distance;
	vector<number> _total_histo, _partial_histo;

public:
	Remoteness();
	virtual ~Remoteness();

	std::string get_output_string(llint curr_step);

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init();
};

extern "C" BaseObservable *make_Remoteness() {
	return new Remoteness();
}

#endif /* REMOTENESS_H_ */
