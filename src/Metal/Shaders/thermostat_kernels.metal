/**
 * @file    thermostat_kernels.metal
 * @brief   Metal kernels for thermostats
 */

#include "shader_utils.h"

using namespace metal;

struct RNGState {
    ulong state;
    ulong inc;
};

// PCG random number generator
uint pcg_next(thread RNGState &rng) {
    ulong oldstate = rng.state;
    rng.state = oldstate * 6364136223846793005UL + (rng.inc | 1);
    uint xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

float pcg_float(thread RNGState &rng) {
    return (float)pcg_next(rng) * 2.3283064e-10f; // / 2^32
}

// Box-Muller transform
void gaussian(thread RNGState &rng, thread m_number &g1, thread m_number &g2) {
    float u1 = pcg_float(rng);
    float u2 = pcg_float(rng);
    
    // Avoid log(0)
    if(u1 < 1e-7) u1 = 1e-7;
    
    float r = sqrt(-2.0f * log(u1));
    float theta = 2.0f * M_PI_F * u2;
    
    g1 = (m_number)(r * cos(theta));
    g2 = (m_number)(r * sin(theta));
}

struct ThermostatArgs {
    m_number rescale_factor;
    m_number pt;
    m_number pr;
    int N;
};

kernel void brownian_thermostat(device RNGState *rng_states [[buffer(0)]],
                                device m_number4 *velocities [[buffer(1)]],
                                device m_number4 *angular_velocities [[buffer(2)]],
                                constant ThermostatArgs &args [[buffer(3)]],
                                uint2 tid [[thread_position_in_grid]]) {
    int idx = tid.x;
    if(idx >= args.N) return;
    
    RNGState rng = rng_states[idx]; // Load state
    
    if(pcg_float(rng) < args.pt) {
        m_number4 v;
        m_number trash;
        m_number vx, vy, vz;
        gaussian(rng, vx, vy);
        gaussian(rng, vz, trash); // Waste one random number or cache it? CUDA code wastes it
        v.x = vx; v.y = vy; v.z = vz;
        
        v.x *= args.rescale_factor;
        v.y *= args.rescale_factor;
        v.z *= args.rescale_factor;
        v.w = 0.0; // Unused or kinetic energy placeholder? CUDA code sets w = KE
        // CUDABrownianThermostat: v.w = (v.x...)*0.5f;
        v.w = (v.x * v.x + v.y * v.y + v.z * v.z) * 0.5f;
        
        velocities[idx] = v;
    }
    
    if(pcg_float(rng) < args.pr) {
        m_number4 L;
        m_number trash;
        m_number Lx, Ly, Lz;
        gaussian(rng, Lx, Ly);
        gaussian(rng, Lz, trash);
        L.x = Lx; L.y = Ly; L.z = Lz;
        
        L.x *= args.rescale_factor;
        L.y *= args.rescale_factor;
        L.z *= args.rescale_factor;
        L.w = (L.x * L.x + L.y * L.y + L.z * L.z) * 0.5f;
        
        angular_velocities[idx] = L;
    }
    
    rng_states[idx] = rng; // Save state
}
