/**
 * @file    AVBDepletion.h
 * @date    30/ott/2018
 * @author  flavio
 *
 */

#ifndef AVBDEPLETION_H_
#define AVBDEPLETION_H_
#include <random>
#include "../../../../../src/Backends/MCMoves/BaseMove.h"

template<typename number>
class AVBDepletion : public BaseMove<number> {
	protected:
		number _delta_rot, _delta_z, _Vb, _rmax, _rmin, _skin, _target_Vb, _histo_dr, _histo_dt, _histo_dz;
		number _sigma_dep, _tryvolume, _mu_dep, _z;
		llint _rejected_n_too_small[2], _rejected_depletion[2], _rejected_overlap[2], _rejected_avb_bias[2], _split_accepted[2];
		int _ntries, _nretries;
		number _delta_max, _tau_max, _zeta_max;

		std::default_random_engine _generator;
		std::poisson_distribution<int> _poisson;

		bool are_close(BaseParticle<number> * , BaseParticle<number> *, bool log);
		number fake_vb (number, number);
		uint32_t _histo[100][100][100], _gain[3];

	public:
		AVBDepletion();
		virtual ~AVBDepletion();

		void apply (llint curr_step);
		virtual void get_settings(input_file &inp, input_file &sim_inp);
		virtual void init(void);
		virtual void log_parameters();
};


extern "C" BaseMove<float> * make_AVBDepletion_float()   { return new AVBDepletion<float> () ; }
extern "C" BaseMove<double> * make_AVBDepletion_double() {return new AVBDepletion<double> () ; }

#endif // AVBDEPLETION_H_
