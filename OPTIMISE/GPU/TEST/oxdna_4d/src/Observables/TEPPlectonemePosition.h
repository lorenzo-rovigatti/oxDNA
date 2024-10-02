/*
 * TEPPlectonemePosition.h
 *
 *  Created on: June 21, 2016
 *      Author: Lorenzo Rovigatti
 */

#ifndef TEPPLECTONEMEPOSITION_H_
#define TEPPLECTONEMEPOSITION_H_

#include "BaseObservable.h"

/**
 * @brief Outputs the position of a plectoneme in a TEP system.
 *
 * To use this observable, use type = TEP_plectoneme_position
 *
 * @verbatim
 [bead_minimum_distance = <int> (the minimum integer separation between beads whose relative distance will be checked. Defaults to 7)]
 [distance_threshold  = <float> (a plectoneme is identified when any two beads with indices farther away than bead_minimum_distance are closer than this distance. Defaults to 2)]
 @endverbatim
 */

class TEPPlectonemePosition: public BaseObservable {
protected:
	int _bead_minimum_distance;
	number _distance_threshold;
	number _old_pos;
	bool _print_pos;

public:
	TEPPlectonemePosition();
	virtual ~TEPPlectonemePosition();

	std::string get_output_string(llint curr_step);
	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
};

#endif /* TEPPLECTONEMEPOSITION_H_ */
