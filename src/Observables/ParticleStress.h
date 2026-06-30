/*
 * ParticleStress.h
 *
 *  Created on: 25/jun/2026
 *      Author: lorenzo
 */
#ifndef PARTICLE_STRESS_H_
#define PARTICLE_STRESS_H_

#include "BaseObservable.h"

/**
 * @brief Outputs the instantaneous particle stress, computed through the virial.
 *
 * The supported syntax is
 * @verbatim
 type = particle_stress (an observable that computes the particle stress of the system)
 @endverbatim
 */
class ParticleStress: public BaseObservable {
protected:
	bool _also_coordinates = false;
public:
	ParticleStress();
	virtual ~ParticleStress();

	void get_settings(input_file &my_inp, input_file &sim_inp);

	std::vector<StressTensor> stress_tensors();
	std::string get_output_string(llint curr_step);
};

#endif /* PARTICLE_STRESS_H_ */
