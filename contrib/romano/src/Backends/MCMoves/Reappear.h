/**
 * @file    Reappear.h
 * @date    20/ott/2018
 * @author  flavio
 *
 *
 */

#ifndef REAPPEAR_H_
#define REAPPEAR_H_
#define MCMOVE_CUSTOM

#include "../../../../../src/Backends/MCMoves/BaseMove.h"

class Reappear : public BaseMove {
	protected:

		LR_vector pos_old;

	public:
		Reappear();
		virtual ~Reappear();

		void apply (llint curr_step);
		virtual void get_settings(input_file &inp, input_file &sim_inp);
		virtual void init(void);
};


extern "C" BaseMove * make_Reappear()   { return new Reappear () ; }

#endif // REAPPEAR_H_
