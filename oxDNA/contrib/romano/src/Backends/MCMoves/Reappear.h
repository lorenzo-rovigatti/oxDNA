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

template<typename number>
class Reappear : public BaseMove<number> {
	protected:

		LR_vector<number> pos_old;

	public:
		Reappear();
		virtual ~Reappear();

		void apply (llint curr_step);
		virtual void get_settings(input_file &inp, input_file &sim_inp);
		virtual void init(void);
};


extern "C" BaseMove<float> * make_Reappear_float()   { return new Reappear<float> () ; }
extern "C" BaseMove<double> * make_Reappear_double() { return new Reappear<double> () ; }

#endif // REAPPEAR_H_
