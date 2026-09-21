/**
 * @file    MetalRNAInteraction.h
 * @brief   Metal RNA/RNA2 interaction adapter
 */

#ifndef METALRNAINTERACTION_H_
#define METALRNAINTERACTION_H_

#include "MetalBaseInteraction.h"

class MetalRNAInteraction: public MetalBaseInteraction {
public:
	MetalRNAInteraction() = default;
	virtual ~MetalRNAInteraction() = default;

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

#endif /* METALRNAINTERACTION_H_ */
