/*
 * ThermostatFactory.cpp
 *
 *  Created on: Feb 15, 2013
 *      Author: Flavio
 */

#include "ThermostatFactory.h"

#include "NoThermostat.h"
#include "BrownianThermostat.h"
#include "BussiThermostat.h"
#include "RefreshThermostat.h"
#include "LangevinThermostat.h"
#include "SRDThermostat.h"
#include "DPDThermostat.h"

#include "../../Utilities/Utils.h"

ThermostatPtr ThermostatFactory::make_thermostat(input_file &inp, BaseBox *box) {
	// The default thermostat is no thermostat
	char thermostat_type[512] = "no";
	getInputString(&inp, "thermostat", thermostat_type, 0);

	if(!strncmp(thermostat_type, "john", 512) || !strncmp(thermostat_type, "brownian", 512)) return std::make_shared<BrownianThermostat>();
	else if(!strncmp(thermostat_type, "no", 512)) return std::make_shared<NoThermostat>();
	else if(!strncmp(thermostat_type, "refresh", 512)) return std::make_shared<RefreshThermostat>();
	else if(!strncmp(thermostat_type, "langevin", 512)) return std::make_shared<LangevinThermostat>();
	else if(!strncmp(thermostat_type, "bussi", 512) || !strncmp(thermostat_type, "Bussi", 512)) return std::make_shared<BussiThermostat>();
	else if(!strncmp(thermostat_type, "srd", 512) || !strncmp(thermostat_type, "SRD", 512)) return std::make_shared<SRDThermostat>(box);
	else if(!strncmp(thermostat_type, "DPD", 512)) return std::make_shared<DPDThermostat>();
	else throw oxDNAException("Invalid Thermostat '%s'", thermostat_type);
}
