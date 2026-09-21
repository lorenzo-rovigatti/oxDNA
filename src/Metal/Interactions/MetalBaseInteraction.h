/**
 * @file    MetalBaseInteraction.h
 * @date    Created for Metal backend
 * @author  antigravity
 */

#ifndef METALBASEINTERACTION_H_
#define METALBASEINTERACTION_H_

#include "../MetalUtils.h"
#include "../Lists/MetalBaseList.h"
#include "../metal_utils/MetalBox.h"
#include "../metal_defs.h"

// Forward declaration
struct Metal_kernel_cfg;

/**
 * @brief Abstract class providing an interface for Metal-based interactions.
 */
class MetalBaseInteraction {
protected:
	bool _use_edge = false;
	bool _edge_compatible = false;
	int _n_forces = 1;
	bool _avoid_cpu_calculations = false;
	bool _use_cpu_fallback = true;
    
    Metal_kernel_cfg _launch_cfg;
    
	bool _update_st = false;
    id<MTLBuffer> _d_st = nil;
    
	int _N = -1;
    
    id<MTLBuffer> _d_edge_forces = nil;
    id<MTLBuffer> _d_edge_torques = nil;
    
    id<MTLDevice> _device;
    id<MTLLibrary> _library;

	virtual void _sum_edge_forces(id<MTLBuffer> d_forces);
	virtual void _sum_edge_forces_torques(id<MTLBuffer> d_forces, id<MTLBuffer> d_torques);

public:
	MetalBaseInteraction();
	virtual ~MetalBaseInteraction();

	virtual void get_settings(input_file &inp) = 0;
	virtual void get_metal_settings(input_file &inp);
	virtual void metal_init(int N, id<MTLDevice> device, id<MTLLibrary> library);
	virtual m_number get_metal_rcut() = 0;

	virtual void compute_forces(MetalBaseList *lists, id<MTLBuffer> d_poss, id<MTLBuffer> d_orientations, id<MTLBuffer> forces,
                       id<MTLBuffer> torques,
                       id<MTLBuffer> bonds,
                       id<MTLBuffer> metal_box,
                       id<MTLBuffer> energies=nil) = 0;

	// Encode the force computation into an existing command buffer without
	// committing or waiting on it. Lets the backend batch a whole MD step into
	// a single GPU submission. The default falls back to the self-contained
	// compute_forces().
	virtual void encode_forces(id<MTLCommandBuffer> command_buffer, MetalBaseList *lists, id<MTLBuffer> d_poss,
	                           id<MTLBuffer> d_orientations, id<MTLBuffer> forces, id<MTLBuffer> torques,
	                           id<MTLBuffer> bonds, id<MTLBuffer> metal_box, id<MTLBuffer> energies = nil) {
		(void) command_buffer;
		compute_forces(lists, d_poss, d_orientations, forces, torques, bonds, metal_box, energies);
	}

	bool use_cpu_fallback() const {
		return _use_cpu_fallback;
	}

	bool avoid_cpu_calculations() const {
		return _avoid_cpu_calculations;
	}
};

#endif /* METALBASEINTERACTION_H_ */
