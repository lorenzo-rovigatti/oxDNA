/**
 * @file    MetalLJInteraction.mm
 * @brief   Metal LJ interaction adapter
 */

#include "MetalLJInteraction.h"

#include "MetalCPUForceFallback.h"
#include "../../Utilities/oxDNAException.h"

void MetalLJInteraction::get_settings(input_file &inp) {
	(void) inp;
}

m_number MetalLJInteraction::get_metal_rcut() {
	return 0.f;
}

void MetalLJInteraction::compute_forces(MetalBaseList *lists,
                                        id<MTLBuffer> d_poss,
                                        id<MTLBuffer> d_orientations,
                                        id<MTLBuffer> d_forces,
                                        id<MTLBuffer> d_torques,
                                        id<MTLBuffer> d_bonds,
                                        id<MTLBuffer> d_box,
                                        id<MTLBuffer> d_energies) {
	(void) lists;
	(void) d_bonds;
	(void) d_box;

	if(!_use_cpu_fallback) {
		throw oxDNAException("Metal LJ native force path is not available. Set Metal_avoid_cpu_calculations = 0");
	}

	MetalCPUForceFallback::compute(_N, d_poss, d_orientations, d_forces, d_torques, d_energies);
}
