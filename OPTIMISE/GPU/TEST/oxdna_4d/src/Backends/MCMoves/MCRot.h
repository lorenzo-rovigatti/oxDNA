/**
 * @file    MCRot.h
 * @date    30/apr/2014
 * @author  flavio
 *
 */

#ifndef MCROT_H_
#define MCROT_H_

#include "BaseMove.h"


class MCRot : public BaseMove {
	protected:
		number _delta;
		LR_matrix _orientation_old, _orientationT_old;

	public:
		MCRot();
		virtual ~MCRot();

		virtual void init();
		virtual void get_settings(input_file &inp, input_file &sim_inp);
		void apply(llint curr_step);
		virtual void log_parameters();
};
#endif
