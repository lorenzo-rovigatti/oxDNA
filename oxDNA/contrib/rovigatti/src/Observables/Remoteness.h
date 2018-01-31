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

template<typename number>
class Remoteness : public BaseObservable<number>  {
protected:
	Cells<number> *_cells;
	BaseParticle<number> *_probe;
	BaseParticle<number> **_particles;
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
	virtual void init(ConfigInfo<number> &config_info);
};

extern "C" BaseObservable<float> *make_float() { return new Remoteness<float>(); }
extern "C" BaseObservable<double> *make_double() { return new Remoteness<double>(); }

#endif /* REMOTENESS_H_ */
