/*
 * ParticleForceAndTorque.h
 *
 *  Created on: May 13, 2025
 *      Author: lorenzo
 */

#ifndef PARTICLEFORCEANDTORQUE_H_
#define PARTICLEFORCEANDTORQUE_H_

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
class ParticleForceAndTorque: public BaseObservable {
protected:
	std::string _particle_string;
	bool _lab_frame = true;
	std::vector<int> _indexes;
	std::set<BaseParticle *> _particles;
public:
	ParticleForceAndTorque();
	virtual ~ParticleForceAndTorque();

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init();
	std::string get_output_string(llint curr_step);
};

#endif /* PARTICLEFORCEANDTORQUE_H_ */
