/*
 * MetalBaseThermostat.mm
 *
 *  Created for Metal backend
 */

#include "MetalBaseThermostat.h"

MetalBaseThermostat::MetalBaseThermostat() {

}

MetalBaseThermostat::~MetalBaseThermostat() {

}

void MetalBaseThermostat::get_settings(input_file &inp) {
    BaseThermostat::get_settings(inp);
}

void MetalBaseThermostat::metal_init(int N, id<MTLDevice> device, id<MTLLibrary> library) {
	_N = N;
    _device = device;
    _library = library;
}
