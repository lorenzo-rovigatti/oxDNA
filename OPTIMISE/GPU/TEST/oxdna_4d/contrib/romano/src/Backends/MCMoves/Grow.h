/**
 * @file    Grow.h
 * @date    20/ott/2018
 * @author  flavio
 *
 *
 */

#ifndef REAPPEAR_H_
#define REAPPEAR_H_
#define MCMOVE_CUSTOM

#include "../../../../../src/Backends/MCMoves/BaseMove.h"
#include "../../../../../src/Particles/SpheroCylinder.h"

template<typename number>
class Grow : public BaseMove<number> {
	protected:
		number _target_sigma;
		std::set<int> _particles_to_grow;
		int _n_to_grow;

	public:
		Grow();
		virtual ~Grow();

		void apply (llint curr_step);
		virtual void get_settings(input_file &inp, input_file &sim_inp);
		virtual void init(void);
		virtual double get_acceptance();
};


extern "C" BaseMove<float> * make_Grow_float()   { return new Grow<float> () ; }
extern "C" BaseMove<double> * make_Grow_double() { return new Grow<double> () ; }

#endif // REAPPEAR_H_
