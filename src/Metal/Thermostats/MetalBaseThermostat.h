/**
 * @file    MetalBaseThermostat.h
 * @date    Created for Metal backend
 * @author  antigravity
 */

#ifndef METALBASETHERMOSTAT_H_
#define METALBASETHERMOSTAT_H_

#include "../MetalUtils.h"
#include "../../Backends/Thermostats/BaseThermostat.h"

// Forward declaration
struct Metal_kernel_cfg;

/**
 * @brief Abstract class providing an interface for Metal-based thermostats.
 */
class MetalBaseThermostat : public virtual BaseThermostat {
protected:
	int _N = -1;
    id<MTLDevice> _device;
    id<MTLLibrary> _library;
    
public:
	MetalBaseThermostat();
	virtual ~MetalBaseThermostat();

	virtual void get_settings(input_file &inp) override;
	virtual void metal_init(int N, id<MTLDevice> device, id<MTLLibrary> library);
    using BaseThermostat::apply; // Unhide overloaded virtual function
    
    // Metal-specific apply
    virtual void apply(id<MTLBuffer> d_velocities, id<MTLBuffer> d_angular_velocities, id<MTLBuffer> d_orientations, id<MTLBuffer> d_forces, id<MTLBuffer> d_torques, id<MTLBuffer> d_poss) = 0;

    // Encode the thermostat into an existing command buffer without committing
    // or waiting. The default falls back to the self-contained apply().
    virtual void encode_apply(id<MTLCommandBuffer> command_buffer, id<MTLBuffer> d_velocities, id<MTLBuffer> d_angular_velocities, id<MTLBuffer> d_orientations, id<MTLBuffer> d_forces, id<MTLBuffer> d_torques, id<MTLBuffer> d_poss) {
        (void) command_buffer;
        apply(d_velocities, d_angular_velocities, d_orientations, d_forces, d_torques, d_poss);
    }
};

#endif /* METALBASETHERMOSTAT_H_ */
