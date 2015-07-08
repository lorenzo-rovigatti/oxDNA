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
 * @brief Outputs either the cosine of the angle between one of the orientation
 * 				vectors of neighbouring particles, or the twist angle in the TEP model, namely
 * 				TODO: write the awful thing
 * 				where u and u' are the longitudinal orientation vectors of neighbouring beads (T.v1 and T.v1')
 * 				and f and f' are the other two orientation vectors.
 *
 * To use this observable, use type = vector_angle
 *
 * This observable takes one mandatory argument and three optional arguments
 * @verbatim
first_particle_index = <int> (defaults to 0. index of the first particle on which to compute the angle with the next particle)
last_particle_index = <int> (defaults to the index of the first-but last bead in the same strand as the first particle. Therefore, if I have a strand of N beads, the last one will be the one with index N-2. This is because the last bead is atypical in the TEP model (e.g. it's aligned with the vector before it rather than the one in front of it.). index of the last particle of the chain on which to compute the angle. )
angle_index = <int> (defaults to 1. Can be 1,2, or 3 depending on the vectors we wish to consider, or 0. In that case it measures the twisting angle, as defined above. Only the twisting angle can be negative, all the others will be always taken positive.
@endverbatim
 */

template<typename number>
class VectorAngle : public BaseObservable<number> {
private:
	// arguments
	int _chain_id;
	int _first_particle_index;
	int _last_particle_index;
	int _angle_index;
	// indices of the first and last particle in the subchain to average
	int _first_particle_id;
	int _last_particle_id;
	//ausiliary variables to compute the cosine - declared here so that
	// they are not declared every time the function get_output_string is called.
public:
	VectorAngle();
	virtual ~VectorAngle();

	virtual void init(ConfigInfo<number> &config_info);
	virtual void get_settings(input_file &my_inp, input_file &sim_inp);

	std::string get_output_string(llint curr_step);
	void set_first_last_particle_id();
};

#endif /* VECTORANGLE_H_ */
