/*
 * PotentialEnergySelected.h
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#ifndef PotentialEnergySelected_H_
#define PotentialEnergySelected_H_

#include "BaseObservable.h"
#include "../Interactions/DNAInteraction.h"

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

class PotentialEnergySelected: public BaseObservable {
protected:
	bool _split;
	std::vector<int> names;
public:
	PotentialEnergySelected();
	virtual ~PotentialEnergySelected();

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);

	std::string get_output_string(llint curr_step);
	number get_potential_energy();
};

#endif /* PotentialEnergySelected_H_ */
