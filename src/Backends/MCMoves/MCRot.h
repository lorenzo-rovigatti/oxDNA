/**
 * @file    MCRot.h
 * @date    30/apr/2014
 * @author  flavio
 *
 */

#ifndef MCROT_H_
#define MCROT_H_

#include "BaseMove.h"

template<typename number>
class MCRot : public BaseMove<number> {
	protected:
		number _delta;
		LR_matrix<number> _orientation_old, _orientationT_old;
	
	public:
		MCRot(ConfigInfo<number> *Info);
		virtual ~MCRot();

		virtual void get_settings(input_file &inp, input_file &sim_inp);
		void apply (llint curr_step);
};
#endif
