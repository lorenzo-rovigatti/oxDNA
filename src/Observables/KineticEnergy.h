/*
 * KineticEnergy.h
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#ifndef KINETICENERGY_H_
#define KINETICENERGY_H_

#include "BaseObservable.h"

#include <set>

/**
 * @brief Outputs the total kinetic energy of the system.
 */

class KineticEnergy: public BaseObservable {
protected:
	std::set<int> _directions;
public:
	KineticEnergy();
	virtual ~KineticEnergy();

	void get_settings(input_file &my_inp, input_file &sim_inp);

	std::string get_output_string(llint curr_step);
	number get_kinetic_energy();
};

#endif /* KINETICENERGY_H_ */
