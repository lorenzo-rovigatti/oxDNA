/**
 * @file    DepletionVolume.h
 * @date    30/ott/2018
 * @author  flavio
 *
 *
 * idea taken from https://arxiv.org/pdf/1508.07077.pdf,
 * but applied rather loosely...
 */

#ifndef DEPLETION_VOLUME_H_
#define DEPLETION_VOLUME_H_
#include <random>
#include "../../../../../src/Backends/MCMoves/BaseMove.h"
#include "../../../../../src/Lists/RodCells.h"

template<typename number>
class DepletionVolume : public BaseMove<number> {
	protected:
		number _delta, _delta_max, _pressure;
		number _sigma_dep, _mu_gas, _z, _length;
	
		bool _change_volume;
		int _affected_side;

		int _ntries;
		std::vector<LR_vector<number> > _pos_old;

		std::default_random_engine _generator;
		std::poisson_distribution<int> _poisson;

	public:
		DepletionVolume();
		virtual ~DepletionVolume();

		void apply (llint curr_step);
		virtual void get_settings(input_file &inp, input_file &sim_inp);
		virtual void init(void);
		virtual void log_parameters();
};


extern "C" BaseMove<float> * make_DepletionVolume_float()   { return new DepletionVolume<float> () ; }
extern "C" BaseMove<double> * make_DepletionVolume_double() { return new DepletionVolume<double> () ; }

#endif // DEPLETION_VOLUME_H_
