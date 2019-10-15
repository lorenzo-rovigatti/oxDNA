/*
 * Pitch.h
 *
 *  Created on: Mar 14, 2014
 *      Author: Ben Snodin
 */

#ifndef PITCH_H_
#define PITCH_H_

#include "BaseObservable.h"

/**
 * @brief Outputs the pitch angle between the specified base-pairs
 *
 * The pitch is defined as the angle between the projection in the plane normal 
 * to the helix axis of the base-base vectors of specified two base pairs
 *
 * @verbatim
 bp1a_id = <int> (base pair 1 particle a id)
 bp1b_id = <int> (base pair 1 particle b id)
 bp2a_id = <int> (base pair 2 particle a id)
 bp2b_id = <int> (base pair 2 particle b id)
 @endverbatim
 */

class Pitch: public BaseObservable {
protected:
	int _bp1a_id;
	int _bp1b_id;
	int _bp2a_id;
	int _bp2b_id;
public:
	Pitch();
	virtual ~Pitch();

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);

	std::string get_output_string(llint curr_step);
};

#endif /* PITCH_H_ */
