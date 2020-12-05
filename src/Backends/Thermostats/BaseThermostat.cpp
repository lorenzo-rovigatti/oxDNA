/*
 * BaseThermostat.cpp
 *
 *  Created on: 14 ott 2019
 *      Author: lorenzo
 */

#include "BaseThermostat.h"

BaseThermostat::BaseThermostat() :
				_T((number) 0.f),
				_supports_shear(false),
				_lees_edwards(false) {
	CONFIG_INFO->subscribe("T_updated", [this]() { this->_on_T_update(); });
}

void BaseThermostat::get_settings(input_file &inp) {
	// we get the temperature from the input file;
	char raw_T[256];
	getInputString(&inp, "T", raw_T, 1);
	_T = Utils::get_temperature(raw_T);

	getInputBool(&inp, "lees_edwards", &_lees_edwards, 0);
	if(_lees_edwards) {
		if(!_supports_shear) throw oxDNAException("The chosen thermostat does not support Lees-Edwards boundary conditions");
		getInputNumber(&inp, "lees_edwards_shear_rate", &_shear_rate, 1);
	}
}

void BaseThermostat::_on_T_update() {
	_T = CONFIG_INFO->temperature();
	init();
}
