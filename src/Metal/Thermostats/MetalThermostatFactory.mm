/*
 * MetalThermostatFactory.mm
 *
 *  Created for Metal backend
 */

#include "MetalThermostatFactory.h"
#include "../../Utilities/oxDNAException.h"

#include "MetalBrownianThermostat.h"
// #include "MetalLangevinThermostat.h" // Pending

MetalBaseThermostat* MetalThermostatFactory::make_thermostat(input_file &inp) {
	std::string thermostat_type;
	if(getInputString(&inp, "thermostat", thermostat_type, 1) == KEY_NOT_FOUND) {
		throw oxDNAException("Cannot search for thermostat in input_file");
	}

	if(thermostat_type.compare("brownian") == 0 || thermostat_type.compare("john") == 0) {
		return new MetalBrownianThermostat();
	}
    // else if(thermostat_type.compare("langevin") == 0) {
    //     return new MetalLangevinThermostat();
    // }
	else if(thermostat_type.compare("no") == 0) {
		return nullptr;
	}
	else {
		throw oxDNAException("Metal thermostat '%s' is not supported yet", thermostat_type.c_str());
	}
}
