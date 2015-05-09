/*
 * ThermostatFactory.cpp
 *
 *  Created on: Feb 15, 2013
 *      Author: Flavio
 */

#include "CUDAThermostatFactory.h"
#include "CUDABrownianThermostat.h"
#include "CUDALangevinThermostat.h"
#include "CUDANoThermostat.h"
#include "CUDASRDThermostat.h"
#include "CUDABussiThermostat.h"
#include "../../Utilities/Utils.h"

CUDAThermostatFactory::CUDAThermostatFactory() {

}

CUDAThermostatFactory::~CUDAThermostatFactory() {

}

template<typename number, typename number4>
CUDABaseThermostat<number, number4> *CUDAThermostatFactory::make_thermostat(input_file &inp, number &box_side) {
	// The default thermostat is no thermostat
	char thermostat_type[512] = "no";
	getInputString(&inp, "thermostat", thermostat_type, 0);

	if(!strncmp(thermostat_type, "john", 512)) return new CUDABrownianThermostat<number, number4>();
	else if(!strncmp(thermostat_type, "brownian", 512)) return new CUDABrownianThermostat<number, number4>();
	else if(!strncmp(thermostat_type, "bussi", 512)) return new CUDABussiThermostat<number, number4>();
	else if(!strncmp(thermostat_type, "langevin", 512)) return new CUDALangevinThermostat<number, number4>();
	else if(!strncmp(thermostat_type, "srd", 512) || !strncmp(thermostat_type, "SRD", 512)) return new CUDASRDThermostat<number, number4>(box_side);
	else if(!strncmp(thermostat_type, "no", 512)) return new CUDANoThermostat<number, number4>();
	else throw oxDNAException("Invalid CUDA thermostat '%s'", thermostat_type);
}

template CUDABaseThermostat<float, float4> *CUDAThermostatFactory::make_thermostat<float, float4>(input_file &inp, float &box_side);
template CUDABaseThermostat<double, LR_double4> *CUDAThermostatFactory::make_thermostat<double, LR_double4>(input_file &inp, double &box_side);
