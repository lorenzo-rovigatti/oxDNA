/*
 * CUDAPatchyInteraction.h
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#ifndef CUDAFSINTERACTION_H_
#define CUDAFSINTERACTION_H_

#include "CUDABaseInteraction.h"

#include "FSInteraction.h"

#define CUDA_MAX_FS_PATCHES 4

/**
 * @brief CUDA implementation of the {@link FSInteraction}.
 */
template<typename number, typename number4>
class CUDAFSInteraction: public CUDABaseInteraction<number, number4>, public FSInteraction<number> {
public:
	CUDAFSInteraction();
	virtual ~CUDAFSInteraction();

	void get_settings(input_file &inp);
	void cuda_init(number box_side, int N);
	number get_cuda_rcut() { return this->get_rcut(); }

	void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds);
};

// CUDA interactions need ad-hoc names, otherwise they will conflict with the entry points for their CPU counterparts
extern "C" IBaseInteraction<float> *make_cuda_float() { return new CUDAFSInteraction<float, float4>(); }
extern "C" IBaseInteraction<double> *make_cuda_double() { return new CUDAFSInteraction<double, LR_double4>(); }

#endif /* CUDAFSINTERACTION_H_ */
