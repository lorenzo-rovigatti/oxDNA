/*
 * MetalBrownianThermostat.mm
 *
 *  Created for Metal backend
 */

#include "MetalBrownianThermostat.h"
#include "../../Utilities/ConfigInfo.h"

// Define RNG state struct to match shader
struct MetalRNGState {
    uint64_t state;
    uint64_t inc;
};

MetalBrownianThermostat::MetalBrownianThermostat() : MetalBaseThermostat(), BrownianThermostat() {
    _d_rng_state = nil;
}

MetalBrownianThermostat::~MetalBrownianThermostat() {
    _d_rng_state = nil;
}

void MetalBrownianThermostat::get_settings(input_file &inp) {
    BrownianThermostat::get_settings(inp);
    MetalBaseThermostat::get_settings(inp);
}

void MetalBrownianThermostat::init() {
    BrownianThermostat::init();
}

void MetalBrownianThermostat::metal_init(int N, id<MTLDevice> device, id<MTLLibrary> library) {
    MetalBaseThermostat::metal_init(N, device, library);
    BrownianThermostat::init();
    
    _init_rng(N);
    
    NSError *error = nil;
    _thermostat_pso = [_device newComputePipelineStateWithFunction:[_library newFunctionWithName:@"brownian_thermostat"] error:&error];
    if(!_thermostat_pso) printf("Error creating brownian_thermostat PSO: %s\n", [[error localizedDescription] UTF8String]);
}

void MetalBrownianThermostat::_init_rng(int N) {
    _d_rng_state = MetalUtils::allocate_buffer<MetalRNGState>(_device, N, MTLResourceStorageModeShared);
    MetalRNGState *states = (MetalRNGState *)[_d_rng_state contents];
    
    // Initialize with simple seed
    srand(1234);
    for(int i=0; i<N; i++) {
        states[i].state = rand();
        states[i].inc = (rand() | 1); // Must be odd
    }
}

void MetalBrownianThermostat::apply(id<MTLBuffer> d_velocities, id<MTLBuffer> d_angular_velocities, id<MTLBuffer> d_orientations, id<MTLBuffer> d_forces, id<MTLBuffer> d_torques, id<MTLBuffer> d_poss) {
    if(!_own_queue) _own_queue = [_device newCommandQueue];
    id<MTLCommandBuffer> commandBuffer = [_own_queue commandBuffer];
    encode_apply(commandBuffer, d_velocities, d_angular_velocities, d_orientations, d_forces, d_torques, d_poss);
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
}

void MetalBrownianThermostat::encode_apply(id<MTLCommandBuffer> commandBuffer, id<MTLBuffer> d_velocities, id<MTLBuffer> d_angular_velocities, id<MTLBuffer> d_orientations, id<MTLBuffer> d_forces, id<MTLBuffer> d_torques, id<MTLBuffer> d_poss) {
    llint step = CONFIG_INFO->curr_step;
    if((step % this->_newtonian_steps) != 0) return;

    id<MTLComputeCommandEncoder> computeEncoder = [commandBuffer computeCommandEncoder];

    [computeEncoder setComputePipelineState:_thermostat_pso];
    [computeEncoder setBuffer:_d_rng_state offset:0 atIndex:0];
    [computeEncoder setBuffer:d_velocities offset:0 atIndex:1];
    [computeEncoder setBuffer:d_angular_velocities offset:0 atIndex:2];
    
    struct {
        m_number rescale_factor;
        m_number pt;
        m_number pr;
        int N;
    } args;
    args.rescale_factor = this->_rescale_factor;
    args.pt = this->_pt;
    args.pr = this->_pr;
    args.N = _N;
    
    [computeEncoder setBytes:&args length:sizeof(args) atIndex:3];
    
    int tpb = 64;
    int blocks = _N / tpb + ((_N % tpb == 0) ? 0 : 1);
    [computeEncoder dispatchThreadgroups:MTLSizeMake(blocks, 1, 1) threadsPerThreadgroup:MTLSizeMake(tpb, 1, 1)];
    
    [computeEncoder endEncoding];
}
