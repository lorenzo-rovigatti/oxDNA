/**
 * @file    MetalBrownianThermostat.h
 * @date    Created for Metal backend
 * @author  antigravity
 */

#ifndef METALBROWNIANTHERMOSTAT_H_
#define METALBROWNIANTHERMOSTAT_H_

#include "MetalBaseThermostat.h"
#include "../../Backends/Thermostats/BrownianThermostat.h"

/**
 * @brief Metal implementation of Brownian thermostat
 */
class MetalBrownianThermostat : public MetalBaseThermostat, public BrownianThermostat {
protected:
    id<MTLBuffer> _d_rng_state;
    id<MTLComputePipelineState> _thermostat_pso;
    id<MTLCommandQueue> _own_queue = nil;

    void _init_rng(int N);

public:
	MetalBrownianThermostat();
	virtual ~MetalBrownianThermostat();

	void get_settings(input_file &inp) override;
    void init() override;
	void metal_init(int N, id<MTLDevice> device, id<MTLLibrary> library) override;
	void apply(id<MTLBuffer> d_velocities, id<MTLBuffer> d_angular_velocities, id<MTLBuffer> d_orientations, id<MTLBuffer> d_forces, id<MTLBuffer> d_torques, id<MTLBuffer> d_poss) override;
	void encode_apply(id<MTLCommandBuffer> command_buffer, id<MTLBuffer> d_velocities, id<MTLBuffer> d_angular_velocities, id<MTLBuffer> d_orientations, id<MTLBuffer> d_forces, id<MTLBuffer> d_torques, id<MTLBuffer> d_poss) override;


    // CPU apply from BaseThermostat - empty for Metal
    void apply(std::vector<BaseParticle *> &particles, llint curr_step) override {}
};

#endif /* METALBROWNIANTHERMOSTAT_H_ */
