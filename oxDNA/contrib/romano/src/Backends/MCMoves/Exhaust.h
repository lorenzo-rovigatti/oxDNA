/**
 * @file    Exhaust.h
 * @date    30/ott/2018
 * @author  flavio
 *
 *
 */

#ifndef EXHAUST_H_
#define EXHAUST_H_

#include "../../../../../src/Backends/MCMoves/BaseMove.h"

template<typename number>
class Exhaust : public BaseMove<number> {
	protected:
		number _delta, _delta_max;
		int _max_clust_size;
		llint _n_found_other;
		llint _n_too_large;

		std::vector<LR_vector<number> > _poss_old;
		std::vector<BaseParticle<number> *> _clust;

	public:
		Exhaust();
		virtual ~Exhaust();

		void apply (llint curr_step);
		virtual void get_settings(input_file &inp, input_file &sim_inp);
		virtual void init(void);
		virtual void log_parameters();
};


extern "C" BaseMove<float> * make_Exhaust_float()   { return new Exhaust<float> () ; }
extern "C" BaseMove<double> * make_Exhaust_double() {return new Exhaust<double> () ; }

#endif // EXHAUST_H_
