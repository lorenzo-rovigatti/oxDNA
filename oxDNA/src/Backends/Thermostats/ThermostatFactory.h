/*
 * ThermostatFactory.h
 *
 *  Created on: Feb 15, 2013
 *      Author: Flavio
 */

#ifndef THERMOSTATFACTORY_H_
#define THERMOSTATFACTORY_H_

#include "BaseThermostat.h"

/**
 * @brief Thermostat factory class.
 *
 * @verbatim
[thermostat = no|refresh|brownian|langevin|srd (Select the simulation thermostat for MD simulations. 'no' means constant-energy simulations. 'refresh' is the Anderson thermostat. 'brownian' is an Anderson-like thermostat that refreshes momenta of randomly chosen particles. 'langevin' implements a regular Langevin thermostat. 'srd' is an (experimental) implementation of a stochastic rotational dynamics algorithm. 'no' and 'brownian' are also available on CUDA. Defaults to 'no'.)]
@endverbatim
 */
class ThermostatFactory {
public:
	ThermostatFactory();
	virtual ~ThermostatFactory();

	/**
	 * @brief Method that returns a pointer to a thermostat object
	 *
	 * This function chooses which object to create and returns it. It
	 * is implemented as a series of if/else clauses
	 *
	 * @param inp input file object reference. The settings of the
	 * thermostat will be taken from it. It usually is an {@link
	 * input_file} object created from an input file
	 * @param box_side simulation box side
	 *
	 * @return the Thermostat Object, which must be of a class derived
	 * from BaseThermostat
	 */
	template<typename number>
	static BaseThermostat<number> *make_thermostat(input_file &inp, number &box_side);
};

#endif /* THERMOSTATFACTORY_H_ */
