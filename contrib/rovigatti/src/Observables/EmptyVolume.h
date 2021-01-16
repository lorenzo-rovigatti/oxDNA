/*
 * EmptyVolume.h
 *
 *  Created on: Jan 14, 2016
 *      Author: Lorenzo Rovigatti
 */

#ifndef EMPTYVOLUME_H_
#define EMPTYVOLUME_H_

#include "Lists/Cells.h"

#include "Observables/BaseObservable.h"
#include "Interactions/DNAInteraction.h"
#include "Utilities/OrderParameters.h"

class EmptyVolume: public BaseObservable {
protected:
	number _probe_diameter;
	number _particle_diameter;
	number _rcut, _sqr_rcut;
	Cells *_cells;
	BaseParticle *_probe;
	std::vector<BaseParticle *> _particles;
	int _N;
	int _tries;

public:
	EmptyVolume();
	virtual ~EmptyVolume();

	std::string get_output_string(llint curr_step);

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init();
};

extern "C" BaseObservable *make_EmptyVolume() {
	return new EmptyVolume();
}

#endif /* CONSTRUCTWISE_H_ */
