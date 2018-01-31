/*
 * ThermostatFactory.cpp
 *
 *  Created on: Feb 15, 2013
 *      Author: Flavio
 */

#include "ThermostatFactory.h"
#include "BrownianThermostat.h"
#include "NoThermostat.h"
#include "RefreshThermostat.h"
#include "LangevinThermostat.h"
#include "SRDThermostat.h"
#include "BussiThermostat.h"
#include "DPDThermostat.h"
#include "../../Utilities/Utils.h"

ThermostatFactory::ThermostatFactory() {

}

ThermostatFactory::~ThermostatFactory() {

}

template<typename number>
BaseThermostat<number> *ThermostatFactory::make_thermostat(input_file &inp, BaseBox<number> * box) {
	// The default thermostat is no thermostat
	char thermostat_type[512] = "no";
	getInputString(&inp, "thermostat", thermostat_type, 0);

	if(!strncmp(thermostat_type, "john", 512) || !strncmp(thermostat_type, "brownian", 512)) return new BrownianThermostat<number>();
	else if(!strncmp(thermostat_type, "no", 512)) return new NoThermostat<number>();
	else if(!strncmp(thermostat_type, "refresh", 512)) return new RefreshThermostat<number>();
	else if(!strncmp(thermostat_type, "langevin", 512)) return new LangevinThermostat<number>();
	else if(!strncmp(thermostat_type, "srd", 512) || !strncmp(thermostat_type, "SRD", 512)) return new SRDThermostat<number>(box);
	else if(!strncmp(thermostat_type, "bussi", 512) || !strncmp(thermostat_type, "Bussi", 512)) return new BussiThermostat<number>();
	else if(!strncmp(thermostat_type, "DPD", 512)) return new DPDThermostat<number>();
	else throw oxDNAException("Invalid Thermostat '%s'", thermostat_type);
}

template BaseThermostat<float> *ThermostatFactory::make_thermostat<float>(input_file &inp, BaseBox<float> * box);
template BaseThermostat<double> *ThermostatFactory::make_thermostat<double>(input_file &inp, BaseBox<double> * box);
