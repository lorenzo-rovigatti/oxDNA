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

template<typename number, typename number4>
class CUDAPolymerInteraction: public CUDABaseInteraction<number, number4>, public PolymerInteraction<number> {
protected:
public:
	CUDAPolymerInteraction();
	virtual ~CUDAPolymerInteraction();

	number get_cuda_rcut() { return this->get_rcut(); }
	void get_settings(input_file &inp);

	void cuda_init(number box_side, int N);

	void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box);
};

extern "C" IBaseInteraction<float> *make_CUDAPolymerInteraction_float() { return new CUDAPolymerInteraction<float, float4>(); }
extern "C" IBaseInteraction<double> *make_CUDAPolymerInteraction_double() { return new CUDAPolymerInteraction<double, LR_double4>(); }

#endif /* CUDAPOLYMERINTERACTION_H_ */
