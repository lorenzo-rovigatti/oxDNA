/*
 * VectorAngle.h
 *
 *  Created on: Dec 9, 2014
 *      Author: Ferdinando Randisi
 */

#ifndef VECTORANGLE_H_
#define VECTORANGLE_H_

#include "BaseObservable.h"
/**
 * @brief Outputs either the cosine of the angle between one of the orientation vectors of neighbouring particles (if angle_index = 1,2,3), OR the amount of twist (measured in turns). 
 *
 * To use this observable, use type = vector_angle.
 *
 * This observable takes four optional arguments
 * @verbatim

 [first_particle_index = <int> (defaults to 0. index of the first particle on which to compute the angle with the next particle.)]
 [last_particle_index = <int> (defaults to the index of the first-but last bead in the same strand as the first particle. Therefore, if I have a strand of N beads, the last one will be the one with index N-2. This is because the last bead is atypical in the TEP model (e.g. it's aligned with the vector before it rather than the one in front of it.). index of the last particle of the chain on which to compute the angle.)]
 [angle_index = <int> (defaults to 1. Can be 1,2, or 3 depending on which orientation vector we want to compute the cosine, or 0. In that case it measures the amount of twist, defined as in the TEP model: (1 + v2.v2'+v3.v3')/(2*pi)*(1 + v1.v1') where v1, v2 and v3 are the orientation vectors of the first particle in a pair and v1', v2' and v3' are the orientation vectors of the following particle.)]
 [print_local_details = <bool> (defaults to true. If true, print the quantity relative to each pair of particles. Otherwise, print their average (for angle_index == 1,2,3) OR their sum (that is, the total twist, for angle_index = 0.)]
 @endverbatim
 */
class VectorAngle: public BaseObservable {
private:
	// arguments
	int _first_particle_index;
	int _last_particle_index;
	int _angle_index;

	bool _print_local_details;
	//ausiliary variables to compute the cosine - declared here so that
	// they are not declared every time the function get_output_string is called.
public:
	VectorAngle();
	virtual ~VectorAngle();

	virtual void init();
	virtual void get_settings(input_file &my_inp, input_file &sim_inp);

	std::string get_output_string(llint curr_step);
};

#endif /* VECTORANGLE_H_ */
