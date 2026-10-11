/**
 * @file    MetalPatchyInteraction.h
 * @brief   Metal patchy interaction adapter
 */

#ifndef METALPATCHYINTERACTION_H_
#define METALPATCHYINTERACTION_H_

#include "MetalBaseInteraction.h"

class MetalPatchyInteraction: public MetalBaseInteraction {
public:
	MetalPatchyInteraction() = default;
	virtual ~MetalPatchyInteraction() = default;

	void get_settings(input_file &inp) override;
	m_number get_metal_rcut() override;
	void compute_forces(MetalBaseList *lists,
	                    id<MTLBuffer> d_poss,
	                    id<MTLBuffer> d_orientations,
	                    id<MTLBuffer> d_forces,
	                    id<MTLBuffer> d_torques,
	                    id<MTLBuffer> d_bonds,
	                    id<MTLBuffer> d_box,
	                    id<MTLBuffer> d_energies) override;
};

#endif /* METALPATCHYINTERACTION_H_ */
