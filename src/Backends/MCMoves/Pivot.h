/*
 * Pivot.h
 *
 *  Created on: Mar 8, 2023
 *      Author: lorenzo
 */

#ifndef SRC_BACKENDS_MCMOVES_PIVOT_H_
#define SRC_BACKENDS_MCMOVES_PIVOT_H_

#include "BaseMove.h"

class Pivot : public BaseMove {
public:
	Pivot();
	virtual ~Pivot();

	void log_parameters();

	void get_settings(input_file &inp, input_file &sim_inp) override;
	void init() override;
	void apply (llint curr_step) override;

protected:
	number _delta = 0.;

	void _copy_particle_details(BaseParticle *dest, BaseParticle *src);
};

#endif /* SRC_BACKENDS_MCMOVES_PIVOT_H_ */
