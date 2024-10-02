/**
 * @file    Swap.h
 * @date    10/lug/2019
 * @author  flavio
 *
 *
 */

#ifndef SWAP_H_
#define SWAP_H_
#define MCMOVE_CUSTOM

#include "../../../../../src/Backends/MCMoves/BaseMove.h"

template<typename number>
class Swap : public BaseMove<number> {
	protected:

	public:
		Swap();
		virtual ~Swap();

		void apply (llint curr_step);
		virtual void get_settings(input_file &inp, input_file &sim_inp);
		virtual void init(void);
};


extern "C" BaseMove<float> * make_Swap_float()   { return new Swap<float> () ; }
extern "C" BaseMove<double> * make_Swap_double() { return new Swap<double> () ; }

#endif // SWAP_H_
