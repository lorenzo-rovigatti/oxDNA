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

class VolumeMove: public BaseMove {
protected:
	number _delta;
	std::vector<LR_vector> _pos_old;

	number _verlet_skin;
	number _P;
	bool _isotropic;

public:
	VolumeMove();
	virtual ~VolumeMove();

	void apply(llint curr_step);
	virtual void init();
	virtual void get_settings(input_file &inp, input_file &sim_inp);
	virtual void log_parameters();
};
#endif // VOLUME_MOVE_H_
