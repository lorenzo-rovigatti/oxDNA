/*
 * PAIRForce.h
 *
 *  Created on: Jun 27, 2014
 *      Author: majid
 */

#ifndef PAIRFORCE_H_
#define PAIRFORCE_H_

#include "BaseObservable.h"

/**
 * @brief Prints out all the non-zero forces in between particles. Select with type=pair_force
 *
 * If the particle_id option is specified, only the forces relative to pair interactions that involve the particle
 * with that id will be printed. To use this observable, use type = pair_force.
 *
 * @verbatim
 [particle_id = <int> (Optional argument. particle id.)]
 @endverbatim
 */

class PairForce: public BaseObservable {
protected:
	int _particle_id;
	bool _print_all_particles;

public:
	PairForce();
	virtual ~PairForce();

	std::string get_output_string(llint curr_step);

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
};

#endif /* PAIRFORCE_H_ */
