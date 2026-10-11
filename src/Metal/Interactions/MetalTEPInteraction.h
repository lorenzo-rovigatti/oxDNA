/**
 * @file    MetalTEPInteraction.h
 * @brief   Metal TEP interaction adapter
 */

#ifndef METALTEPINTERACTION_H_
#define METALTEPINTERACTION_H_

#include "MetalBaseInteraction.h"

class MetalTEPInteraction: public MetalBaseInteraction {
public:
	MetalTEPInteraction() = default;
	virtual ~MetalTEPInteraction() = default;

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

#endif /* METALTEPINTERACTION_H_ */
