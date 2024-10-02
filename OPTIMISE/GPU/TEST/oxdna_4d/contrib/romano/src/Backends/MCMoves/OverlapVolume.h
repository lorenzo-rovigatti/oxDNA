/**
 * @file    OverlapVolume.h
 * @date    30/apr/2014
 * @author  flavio
 *
 *
 */

#ifndef MCMOVEOVERLAPVOLUME_H_
#define MCMOVEOVERLAPVOLUME_H_
#define MCMOVE_CUSTOM

#include "../../../../../src/Backends/MCMoves/BaseMove.h"

template<typename number>
class OverlapVolume : public BaseMove<number> {
	protected:

		llint _print_every;
		number _rmax;

		llint _noverlaps;

		number _mesh_r, _mesh_w;
		std::vector<std::vector<llint> > _suc, _try;

		LR_vector<number> _pos_old;
		LR_matrix<number> _orientation_old;
		LR_matrix<number> _orientationT_old;

	public:
		OverlapVolume();
		virtual ~OverlapVolume();

		void apply (llint curr_step);
		virtual void get_settings(input_file &inp, input_file &sim_inp);
		virtual void init(void);

};


extern "C" BaseMove<float> * make_OverlapVolume_float()   { return new OverlapVolume<float> () ; }
extern "C" BaseMove<double> * make_OverlapVolume_double() {return new OverlapVolume<double> () ; }

#endif // MCMOVEOVERLAPVOLUME_H_
