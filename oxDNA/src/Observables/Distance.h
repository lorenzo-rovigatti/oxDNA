/*
 * Distance.h
 *
 *  Created on: Jan 22, 2014
 *      Author: Flavio
 */

#ifndef DISTANCE_H_
#define DISTANCE_H_

#include "BaseObservable.h"

/**
 * @brief Outputs the distance between two particles; can be projected along a direction.
 *
 * To use this observable, use type = distance 
 *
 * This observable takes one mandatory arguments and an tow optional ones:
 * @verbatim
particle_1 = <int> (index of the first particle)
particle_2 = <int> (index of the second particle. The distance is returned as r(2) - r(1))
[PBC = <bool> (Whether to honour PBC. Defaults to True)]
[dir = <float>, <float>, <float> (vector to project the distance along. Beware that it gets normalized after reading. Defaults to (1, 1, 1) / sqrt(3))]
@endverbatim
 */
template<typename number>
class Distance : public BaseObservable<number> {
private:
	int _i, _j;
        bool _PBC, _have_dir;
	LR_vector<number> _dir;
public:
	Distance();
	virtual ~Distance();

	virtual void init(ConfigInfo<number> &config_info); 
	virtual std::string get_output_string(llint curr_step);
	void get_settings (input_file &my_inp, input_file &sim_inp);
};

#endif /* STEP_H_ */
