/*
 * CUDAFSInteraction.h
 *
 *  Created on: 10/apr/2015
 *      Author: lorenzo
 */

#ifndef CUDAFSINTERACTION_H_
#define CUDAFSINTERACTION_H_

#include "CUDA/Interactions/CUDABaseInteraction.h"

#include "FSInteraction.h"

struct CUDA_FS_bonding_pattern;
struct swap_event;

/**
 * @brief CUDA implementation of the {@link FSInteraction}.
 */
template<typename number, typename number4>
class CUDAFSInteraction: public CUDABaseInteraction<number, number4>, public FSInteraction<number> {
protected:
	bool _analyse_bonds;
	bool _called_once;
	CUDA_FS_bonding_pattern *_d_bonds, *_d_new_bonds;
	swap_event *_h_events, *_d_events;
	number4 *_d_three_body_forces, *_d_three_body_torques;
	llint _step;

	// polymer stuff
	int *_d_bonded_neighs;
public:
	CUDAFSInteraction();
	virtual ~CUDAFSInteraction();

	void get_settings(input_file &inp);
	void cuda_init(number box_side, int N);
	number get_cuda_rcut() { return this->get_rcut(); }

	void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box);
};

extern "C" IBaseInteraction<float> *make_CUDAFSInteraction_float() { return new CUDAFSInteraction<float, float4>(); }
extern "C" IBaseInteraction<double> *make_CUDAFSInteraction_double() { return new CUDAFSInteraction<double, LR_double4>(); }

#endif /* CUDAFSINTERACTION_H_ */
