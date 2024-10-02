/*
 * ParticlePosition.h
 *
 *  Created on: Feb 12, 2013
 *      Author: petr
 */

#ifndef PARTICLEPOSITION_H_
#define PARTICLEPOSITION_H_

#include "BaseObservable.h"

/**
 * @brief Outputs the position of the center of mass of a particle.
 *
 * @verbatim
 particle_id = <int> (particle id)
 [orientation = <bool> (defaults to false. If 1, it also prints out the orientation)]
 [absolute = <bool> (defaults to false. If 1, does not use periodic boundaries and it prints out the absolute position of the center of mass)]
 @endverbatim
 */
class ParticlePosition: public BaseObservable {
protected:
	int _particle_id;
	bool _orientation;
	bool _absolute;
public:
	ParticlePosition();
	virtual ~ParticlePosition();

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init();
	std::string get_output_string(llint curr_step);
	//number get_potential_energy();
};

#endif /* POTENTIALENERGY_H_ */
