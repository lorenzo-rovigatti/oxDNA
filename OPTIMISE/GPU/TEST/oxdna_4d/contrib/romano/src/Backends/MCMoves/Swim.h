/**
 * @file    Swim.h
 * @date    30/apr/2014
 * @author  flavio
 *
 *
 */

#ifndef MCMOVEPS_H_
#define MCMOVEPS_H_
#define MCMOVE_CUSTOM

#include "../../../../../src/Backends/MCMoves/BaseMove.h"

template<typename number>
class Swim : public BaseMove<number> {
	protected:
		number _delta;
		int _axis_i;

		LR_vector<number> * _get_axis(BaseParticle<number> * p);

		LR_vector<number> pos_old;

	public:
		Swim();
		virtual ~Swim();

		// we can use any of the three orientation vectors
		enum {
			O_V1 = 0,
			O_V2 = 1,
			O_V3 = 2
		};


		void apply (llint curr_step);
		virtual void get_settings(input_file &inp, input_file &sim_inp);
		virtual void init(void);
		virtual void log_parameters();
};


extern "C" BaseMove<float> * make_Swim_float()   { return new Swim<float> () ; }
extern "C" BaseMove<double> * make_Swim_double() {return new Swim<double> () ; }

#endif // Swim.h
