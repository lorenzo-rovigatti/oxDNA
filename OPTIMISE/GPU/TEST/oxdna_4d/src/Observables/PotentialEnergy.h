/*
 * PotentialEnergy.h
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#ifndef POTENTIALENERGY_H_
#define POTENTIALENERGY_H_

#include "BaseObservable.h"

/**
 * @brief Outputs the potential energy of the system.
 *
 * For now, it outputs the total potential energy or the different contributions to the potential energy.
 * All the energies are per particle.
 *
 * @verbatim
 [split = <bool> (defaults to false, it tells the observable to print all the terms contributing to the potential energy)]
 @endverbatim
 */

class PotentialEnergy: public BaseObservable {
protected:
	bool _split;
public:
	PotentialEnergy();
	virtual ~PotentialEnergy();

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);

	std::string get_output_string(llint curr_step);
	number get_potential_energy();
};

#endif /* POTENTIALENERGY_H_ */
