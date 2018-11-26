/**
 * @file    NDepletion.h
 * @date    30/ott/2018
 * @author  flavio
 *
 *
 * idea taken from https://arxiv.org/pdf/1508.07077.pdf,
 * but applied rather loosely...
 */

#ifndef DEPLETION_H_
#define DEPLETION_H_
#include <random>
#include "../../../../../src/Backends/MCMoves/BaseMove.h"

template<typename number>
class NDepletion : public BaseMove<number> {
	protected:
		number _sigma_dep, _tryvolume, _mu_dep, _z;
		int _ntries;

		std::default_random_engine _generator;
		std::poisson_distribution<int> _poisson;

	public:
		NDepletion();
		virtual ~NDepletion();

		void apply (llint curr_step);
		virtual void get_settings(input_file &inp, input_file &sim_inp);
		virtual void init(void);
		virtual void log_parameters();
};


extern "C" BaseMove<float> * make_NDepletion_float()   { return new NDepletion<float> () ; }
extern "C" BaseMove<double> * make_NDepletion_double() {return new NDepletion<double> () ; }

#endif // DEPLETION_H_
