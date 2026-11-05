/**
 * @file    MetalTEPInteraction.mm
 * @brief   Metal TEP interaction adapter
 */

#include "MetalTEPInteraction.h"

#include "MetalCPUForceFallback.h"
#include "../../Utilities/oxDNAException.h"

void MetalTEPInteraction::get_settings(input_file &inp) {
	(void) inp;
}

m_number MetalTEPInteraction::get_metal_rcut() {
	return 0.f;
}

void MetalTEPInteraction::compute_forces(MetalBaseList *lists,
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
		throw oxDNAException("Metal TEP native force path is not available. Set Metal_avoid_cpu_calculations = 0");
	}

	MetalCPUForceFallback::compute(_N, d_poss, d_orientations, d_forces, d_torques, d_energies);
}
