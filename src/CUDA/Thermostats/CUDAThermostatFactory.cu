/*
 * ThermostatFactory.cpp
 *
 *  Created on: Feb 15, 2013
 *      Author: Flavio
 */

#include "CUDAThermostatFactory.h"

#include "CUDANoThermostat.h"
/*#include "CUDABrownianThermostat.h"
 #include "CUDALangevinThermostat.h"
 #include "CUDASRDThermostat.h"
 #include "CUDABussiThermostat.h"*/
#include "../../Utilities/Utils.h"

CUDABaseThermostat*CUDAThermostatFactory::make_thermostat(input_file &inp, BaseBox * box) {
	// The default thermostat is no thermostat
	char thermostat_type[512] = "no";
	getInputString(&inp, "thermostat", thermostat_type, 0);

	/*if(!strncmp(thermostat_type, "john", 512)) return new CUDABrownianThermostat();
	 else if(!strncmp(thermostat_type, "brownian", 512)) return new CUDABrownianThermostat();
	 else if(!strncmp(thermostat_type, "bussi", 512)) return new CUDABussiThermostat();
	 else if(!strncmp(thermostat_type, "langevin", 512)) return new CUDALangevinThermostat();
	 else if(!strncmp(thermostat_type, "srd", 512) || !strncmp(thermostat_type, "SRD", 512)) {
	 return new CUDASRDThermostat(box);
	 }
	 else */if(!strncmp(thermostat_type, "no", 512)) return new CUDANoThermostat();
	else throw oxDNAException("Invalid CUDA thermostat '%s'", thermostat_type);
}
