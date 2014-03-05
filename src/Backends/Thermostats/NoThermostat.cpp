/*
 * NoThermostat.cpp
 *
 *  Created on: Feb 15, 2013
 *      Author: Flavio
 */

#include "NoThermostat.h"

template<typename number>
NoThermostat<number>::NoThermostat () : BaseThermostat<number>(){
}

template<typename number>
NoThermostat<number>::~NoThermostat () {

}


template<typename number>
void NoThermostat<number>::apply (BaseParticle<number> **particles, llint curr_step) {
	return;
}

template class NoThermostat<float>;
template class NoThermostat<double>;

