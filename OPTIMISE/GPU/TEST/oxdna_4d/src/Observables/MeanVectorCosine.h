/*
 * MeanVectorCosine.h
 *
 *  Created on: Dec 9, 2014
 *      Author: Ferdinando Randisi
 */

#ifndef MEANVECTORCOSINE_H_
#define MEANVECTORCOSINE_H_

#include "BaseObservable.h"
/**
 * @brief Outputs either the mean cosine of the angle between one of the orientation
 * 				vectors of neighbouring particles, or the quantity
 * 				(v2*v2')(v3*v3') - |v2 ^ v2||v3 ^ v3|
 * 				where vi and vi' are the orientation vectors of neighbouring particles.
 *
 * To use this observable, use type = mean_vector_cosine
 *
 * This observable takes one mandatory argument and three optional arguments
 * @verbatim
 chain_id = <int> (chain id)
 first_particle_position = <int> (defaults to 0. position along the chain of  the first particle on which to compute the vector's cosine with the next particle)
 last_particle_position = <int> (defaults to N-2, where N is the number of elements of the chain. Position along the chain of the last particle over which to compute the vector's cosine with the next particle)
 vector_to_average = <int> (defaults to 1. Can be 1,2, or 3 depending on the vectors we wish to consider, or 0. In that case it measures the quantity (v2*v2')(v3*v3') - |v2 ^ v2||v3 ^ v3|)
 @endverbatim
 */

class MeanVectorCosine: public BaseObservable {
private:
	// arguments
	int _chain_id;
	int _first_particle_position;
	int _last_particle_position;
	int _vector_to_average;
	// indices of the first and last particle in the subchain to average
	int _first_particle_id;
	int _last_particle_id;
	number _one_over_number_of_values;
	//ausiliary variables to compute the cosine - declared here so that
	// they are not declared every time the function get_output_string is called.
	LR_vector _u, _up, _v, _vp, _f, _fp;
public:
	MeanVectorCosine();
	virtual ~MeanVectorCosine();

	virtual void init();
	virtual void get_settings(input_file &my_inp, input_file &sim_inp);

	std::string get_output_string(llint curr_step);
	void set_first_last_particle_id();
};

#endif /* MEANVECTORCOSINE_H_ */
