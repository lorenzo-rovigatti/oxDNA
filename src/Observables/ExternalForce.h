/*
 * ExternalForce.h
 *
 *  Created on: Feb 02, 2024
 *      Author: rovigatti
 */

#ifndef EXTERNALFORCE_H_
#define EXTERNALFORCE_H_

#include "BaseObservable.h"

/**
 * @brief Outputs the force applied by one or more external forces.
 */

class ExternalForce: public BaseObservable {
protected:
	std::string _particle_ids;
	std::vector<BaseParticle *> _particles;
	std::vector<LR_vector> _force_averages;
public:
	ExternalForce();
	virtual ~ExternalForce();

	void get_settings(input_file &my_inp, input_file &sim_inp) override;
	void init() override;

	void update_data(llint curr_step) override;

	std::string get_output_string(llint curr_step);
};

#endif /* EXTERNALFORCE_H_ */
