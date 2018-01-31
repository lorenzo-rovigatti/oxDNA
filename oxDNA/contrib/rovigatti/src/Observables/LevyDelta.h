/*
 * LevyDelta.h
 *
 *  Created on: Jan 14, 2016
 *      Author: Lorenzo Rovigatti
 */

#ifndef LEVY_H_
#define LEVY_H_

#include "Lists/Cells.h"

#include "Observables/BaseObservable.h"
#include "Interactions/DNAInteraction.h"

#include "../Interactions/LevyInteraction.h"

template<typename number>
class LevyDelta : public BaseObservable<number>  {
protected:
	vector<BaseParticle<number> *> _tetramer;
	vector<BaseParticle<number> *> _dimer;
	int _tries;
	number _temperature;
	number _patchy_rcut;
	LevyInteraction<number> *_inter;

	BaseParticle<number> *_get_random_particle(vector<BaseParticle<number> *> &, int);
	void _rototranslate(vector<BaseParticle<number> *> &, LR_vector<number> &, BaseParticle<number> *, LR_matrix<number> &);

public:
	LevyDelta();
	virtual ~LevyDelta();

	std::string get_output_string(llint curr_step);

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init(ConfigInfo<number> &config_info);
};

extern "C" BaseObservable<float> *make_float() { return new LevyDelta<float>(); }
extern "C" BaseObservable<double> *make_double() { return new LevyDelta<double>(); }

#endif /* LEVY_H_ */
