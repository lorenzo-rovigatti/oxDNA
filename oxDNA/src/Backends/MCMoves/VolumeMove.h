/**
 * @file    Volume.h
 * @date    30/may/2015
 * @author  flavio
 *
 *
 */

#ifndef VOLUME_MOVE_H_
#define VOLUME_MOVE_H_

#include "BaseMove.h"

template<typename number>
class VolumeMove : public BaseMove<number> {
	protected:
		number _delta;
		std::vector<LR_vector<number> > _pos_old;

		number _verlet_skin;
		number _P;
	
	public:
		VolumeMove(ConfigInfo<number> *Info);
		virtual ~VolumeMove();

		void apply (llint curr_step);
		virtual void init ();
		virtual void get_settings(input_file &inp, input_file &sim_inp);
};
#endif // VOLUME_MOVE_H_
