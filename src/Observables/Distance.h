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
 * @brief Outputs the distance between two particles or centres of mass of sets of particles; can be projected along a direction.
 *
 * To use this observable, use type = distance 
 *
 * This observable takes two mandatory arguments and two optional ones:
 * @verbatim
 particle_1 = <int> (index of the first particle or comma-separated list of particle indexes composing the first set)
 particle_2 = <int> (index of the second particle or comma-separated list of the particle indexes composing the second set. The distance is returned as r(2) - r(1))
 [PBC = <bool> (Whether to honour PBC. Defaults to True)]
 [dir = <float>, <float>, <float> (vector to project the distance along. Beware that it gets normalized after reading. Defaults to (1, 1, 1) / sqrt(3))]
 @endverbatim
 */

class Distance: public BaseObservable {
private:
	std::string _p1_string;
	std::string _p2_string;

	std::set<BaseParticle *> _p1_list;
	std::set<BaseParticle *> _p2_list;

	bool _PBC, _have_dir;
	LR_vector _dir;

	void _check_index(int idx, int N) {
		if(idx < 0 || idx >= N) throw oxDNAException("Distance: invalid id %d", idx);
	}
public:
	Distance();
	virtual ~Distance();

	virtual void init();
	virtual std::string get_output_string(llint curr_step);
	void get_settings(input_file &my_inp, input_file &sim_inp);
};

#endif /* STEP_H_ */
