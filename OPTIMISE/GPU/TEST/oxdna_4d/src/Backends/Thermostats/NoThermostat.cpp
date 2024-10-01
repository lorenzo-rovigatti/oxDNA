/*
 * NoThermostat.cpp
 *
 *  Created on: Feb 15, 2013
 *      Author: Flavio
 */

#include "NoThermostat.h"

NoThermostat::NoThermostat() :
				BaseThermostat() {
}

NoThermostat::~NoThermostat() {

}

void NoThermostat::apply(std::vector<BaseParticle*> &particles, llint curr_step) {
	return;
}

