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

template<typename number>
class ShapeMove : public BaseMove<number> {
	protected:
		number _delta;
		std::vector<LR_vector<number> > _pos_old;

		number _verlet_skin;
		number _P;
	
	public:
		ShapeMove(ConfigInfo<number> *Info);
		virtual ~ShapeMove();

		void apply (llint curr_step);
		virtual void init ();
		virtual void get_settings(input_file &inp, input_file &sim_inp);
};
#endif // SHAPE_MOVE_H_
