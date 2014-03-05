/*
 * KineticEnergy.h
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#ifndef KINETICENERGY_H_
#define KINETICENERGY_H_

#include "BaseObservable.h"

/**
 * @brief Outputs the total kinetic energy of the system.
 */
template<typename number>
class KineticEnergy : public BaseObservable<number> {
public:
	KineticEnergy();
	virtual ~KineticEnergy();

	std::string get_output_string(llint curr_step);
	number get_kinetic_energy();
};

#endif /* KINETICENERGY_H_ */
