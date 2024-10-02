/**
 * @file    ShapeMove.h
 * @date    17/nov/2015
 * @author  flavio
 *
 *
 */

#ifndef SHAPE_MOVE_H_
#define SHAPE_MOVE_H_

#include "BaseMove.h"


class ShapeMove : public BaseMove {
	protected:
		number _delta;
		std::vector<LR_vector > _pos_old;

		number _verlet_skin;
		number _P;

	public:
		ShapeMove();
		virtual ~ShapeMove();

		void apply (llint curr_step);
		virtual void init ();
		virtual void get_settings(input_file &inp, input_file &sim_inp);
		virtual void log_parameters();
};
#endif // SHAPE_MOVE_H_
