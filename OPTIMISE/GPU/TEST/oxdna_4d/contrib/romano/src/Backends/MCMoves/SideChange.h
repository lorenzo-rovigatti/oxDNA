/**
 * @file    SideChange.h
 * @date    25/ott/2018
 * @author  flavio
 *
 *
 */

#ifndef SIDECHANGE_H_
#define SIDECHANGE_H_

#include "../../../../../src/Backends/MCMoves/BaseMove.h"

class SideChange : public BaseMove {
	protected:
		bool _apply_x;
		bool _apply_y;
		bool _apply_z;

		number _delta;
		std::vector<LR_vector> _pos_old;

		number _verlet_skin;
		number _P;

	public:
		SideChange();
		virtual ~SideChange();

		void apply (llint curr_step);
		virtual void init ();
		virtual void get_settings(input_file &inp, input_file &sim_inp);
		virtual void log_parameters();
};

extern "C" BaseMove * make_SideChange()   { return new SideChange () ; }

#endif // SIDECHANGE_H_
