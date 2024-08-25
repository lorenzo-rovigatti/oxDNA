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
	std::vector<BaseParticle *> _tetramer;
	std::vector<BaseParticle *> _dimer;
	int _tries;
	number _temperature;
	number _patchy_rcut;
	LevyInteraction *_inter;

	BaseParticle *_get_random_particle(std::vector<BaseParticle *> &, int);
	void _rototranslate(std::vector<BaseParticle *> &, LR_vector &, BaseParticle *, LR_matrix &);

public:
	LevyDelta();
	virtual ~LevyDelta();

	std::string get_output_string(llint curr_step);

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init();
};

extern "C" BaseObservable *make_LevyDelta() {
	return new LevyDelta();
}

#endif /* LEVY_H_ */
