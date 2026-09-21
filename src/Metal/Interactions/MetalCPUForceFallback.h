/**
 * @file    MetalCPUForceFallback.h
 * @brief   CPU force fallback helper for Metal interactions
 */

#ifndef METALCPUFORCEFALLBACK_H_
#define METALCPUFORCEFALLBACK_H_

#include "../MetalUtils.h"

class MetalCPUForceFallback {
public:
	static void compute(int N,
	                    id<MTLBuffer> d_poss,
	                    id<MTLBuffer> d_orientations,
	                    id<MTLBuffer> d_forces,
	                    id<MTLBuffer> d_torques,
	                    id<MTLBuffer> d_energies = nil);
};

#endif /* METALCPUFORCEFALLBACK_H_ */
