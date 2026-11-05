/*
 * MetalThermostatFactory.h
 *
 *  Created for Metal backend
 */

#ifndef METALTHERMOSTATFACTORY_H_
#define METALTHERMOSTATFACTORY_H_

#include "MetalBaseThermostat.h"

/**
 * @brief Static factory class for Metal thermostats
 */
class MetalThermostatFactory {
public:
	MetalThermostatFactory() = delete;
	virtual ~MetalThermostatFactory() = delete;

	static MetalBaseThermostat *make_thermostat(input_file &inp);
};

#endif /* METALTHERMOSTATFACTORY_H_ */
