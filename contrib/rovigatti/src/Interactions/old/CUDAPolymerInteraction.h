/*
 * CUDAPolymerInteraction.h
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#ifndef CUDAPOLYMERINTERACTION_H_
#define CUDAPOLYMERINTERACTION_H_

#include "CUDA/Interactions/CUDABaseInteraction.h"
#include "PolymerInteraction.h"

class CUDAPolymerInteraction: public CUDABaseInteraction, public PolymerInteraction {
protected:
public:
	CUDAPolymerInteraction();
	virtual ~CUDAPolymerInteraction();

	c_number get_cuda_rcut() {
		return this->get_rcut();
	}
	void get_settings(input_file &inp) override;

	void cuda_init(int N) override;

	void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box);
};

extern "C" BaseInteraction *make_CUDAPolymerInteraction() {
	return new CUDAPolymerInteraction();
}

#endif /* CUDAPOLYMERINTERACTION_H_ */
