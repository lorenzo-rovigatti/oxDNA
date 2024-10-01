/*
 * Density.h
 *
 *  Created on: Oct 30, 2013
 *      Author: Flavio
 */

#ifndef DENSITY_H_
#define DENSITY_H_

#include "BaseObservable.h"

/**
 * @brief Outputs the density of the system in number of particles over volume in simulation units
 *
 */

class Density: public BaseObservable {
public:
	Density();
	virtual ~Density();

	virtual std::string get_output_string(llint curr_step);
};

#endif /* DENSITY_H_ */
