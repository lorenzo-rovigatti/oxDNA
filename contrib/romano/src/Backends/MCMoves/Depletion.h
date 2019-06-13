/**
 * @file    Depletion.h
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

#define DEPLETION_TRANSLATION (0)
#define DEPLETION_ROTATION (1)
#define DEPLETION_SWIM (2)

template<typename number>
class Depletion : public BaseMove<number> {
	protected:
		number _delta_trs, _delta_trs_max, _delta_rot, _delta_rot_max, _delta_swm, _delta_swm_max;
		number _length, _sigma_dep, _rho_dep, _tryvolume, _mu_gas, _z;
		int _ntries, _avg_dn, _n_dn, _n_dn_n0;

		std::default_random_engine _generator;
		std::poisson_distribution<int> _poisson;

	public:
		Depletion();
		virtual ~Depletion();

		void apply (llint curr_step);
		virtual void get_settings(input_file &inp, input_file &sim_inp);
		virtual void init(void);
		virtual void log_parameters();
};


extern "C" BaseMove<float> * make_Depletion_float()   { return new Depletion<float> () ; }
extern "C" BaseMove<double> * make_Depletion_double() {return new Depletion<double> () ; }

#endif // DEPLETION_H_
