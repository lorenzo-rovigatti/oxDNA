/**
 * @file    CutVolume.h
 * @date    25/ott/2018
 * @author  flavio
 *
 *
 */

#ifndef CUTVOLUME_MOVE_H_
#define CUTVOLUME_MOVE_H_

#include "../../../../../src/Backends/MCMoves/BaseMove.h"

template<typename number>
class CutVolume : public BaseMove<number> {
	protected:
		number _delta;
		std::vector<LR_vector<number> > _pos_old;

		number _verlet_skin;
		number _P;

	public:
		CutVolume();
		virtual ~CutVolume();

		void apply (llint curr_step);
		virtual void init ();
		virtual void get_settings(input_file &inp, input_file &sim_inp);
		virtual void log_parameters();
};

extern "C" BaseMove<float> * make_CutVolume_float()   { return new CutVolume<float> () ; }
extern "C" BaseMove<double> * make_CutVolume_double() { return new CutVolume<double> () ; }

#endif // CUTVOLUME_MOVE_H_
