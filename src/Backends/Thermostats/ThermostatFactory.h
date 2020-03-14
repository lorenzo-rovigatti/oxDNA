/*
 * ThermostatFactory.h
 *
 *  Created on: Feb 15, 2013
 *      Author: Flavio
 */

#ifndef THERMOSTATFACTORY_H_
#define THERMOSTATFACTORY_H_

#include "BaseThermostat.h"
#include "../../Boxes/BaseBox.h"

/**
 * @brief Thermostat factory class.
 *
 * @verbatim
[thermostat = no|refresh|brownian|langevin|srd|bussi (Select the simulation thermostat for MD simulations. 'no' means constant-energy simulations. 'refresh' is the Andersen thermostat. 'brownian' is an Anderson-like thermostat that refreshes momenta of randomly chosen particles. 'langevin' implements a regular Langevin thermostat. 'srd' is an (experimental) implementation of a stochastic rotational dynamics algorithm, 'bussi' is the Bussi-Donadio-Parrinello thermostat. 'no', 'brownian' and 'bussi' are also available on CUDA. Defaults to 'no'.)]
@endverbatim
 */
class ThermostatFactory {
public:
	ThermostatFactory() = delete;
	virtual ~ThermostatFactory() = delete;

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
	
	static ThermostatPtr make_thermostat(input_file &inp, BaseBox *box);
};

#endif /* THERMOSTATFACTORY_H_ */
