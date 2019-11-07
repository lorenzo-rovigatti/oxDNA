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

class LevyDelta: public BaseObservable {
protected:
	vector<BaseParticle *> _tetramer;
	vector<BaseParticle *> _dimer;
	int _tries;
	number _temperature;
	number _patchy_rcut;
	LevyInteraction *_inter;

	BaseParticle *_get_random_particle(vector<BaseParticle *> &, int);
	void _rototranslate(vector<BaseParticle *> &, LR_vector &, BaseParticle *, LR_matrix &);

public:
	LevyDelta();
	virtual ~LevyDelta();

	std::string get_output_string(llint curr_step);

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init(ConfigInfo &config_info);
};

extern "C" BaseObservable *make_LevyDelta() {
	return new LevyDelta();
}

#endif /* LEVY_H_ */
