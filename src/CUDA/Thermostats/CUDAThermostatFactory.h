/*
 * CUDAThermostatFactory.h
 *
 *  Created on: Feb 15, 2013
 *      Author: Flavio
 */

#ifndef CUDATHERMOSTATFACTORY_H_
#define CUDATHERMOSTATFACTORY_H_

#include "CUDABaseThermostat.h"

/**
 * @brief CUDA Thermostat factory class.
 */
class CUDAThermostatFactory {
public:
	CUDAThermostatFactory();
	virtual ~CUDAThermostatFactory();

	/**
	 * @brief Method that returns a pointer to a CUDABaseThermostat object
	 *
	 * This function chooses which object to create and returns it. It
	 * is implemented as a series of if/else clauses
	 *
	 * @param inp input file object reference. The settings of the
	 * thermostat will be taken from it. It usually is an input_file object created from an input file
	 * @param box_side
	 *
	 * @return a pointer to a CUDA thermostat Object, which must be of a class derived
	 * from BaseThermostat
	 */
	template<typename number, typename number4>
	static CUDABaseThermostat<number, number4> *make_thermostat(input_file &inp, number &box_side);
};

#endif /* CUDATHERMOSTATFACTORY_H_ */
