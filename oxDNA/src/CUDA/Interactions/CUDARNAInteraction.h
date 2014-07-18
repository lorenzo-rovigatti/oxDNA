/*
 * CUDARNAInteraction.h
 *
 *  Created on: 22/jul/2014
 *      Author: petr
 */

#ifndef CUDARNAINTERACTION_H_
#define CUDARNAINTERACTION_H_

#include "CUDABaseInteraction.h"
#include "../../Interactions/RNAInteraction.h"

/**
 * @brief CUDA implementation of the oxRNA model, as provided by RNAInteraction.
 */
template<typename number, typename number4>
class CUDARNAInteraction: public CUDABaseInteraction<number, number4>, public RNAInteraction<number> {
public:

	bool _grooving;
	CUDARNAInteraction();
	virtual ~CUDARNAInteraction();

	void get_settings(input_file &inp);
	void cuda_init(number box_side, int N);
	number get_cuda_rcut() { return this->get_rcut(); }

	void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, LR_GPU_matrix<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds);
};

#endif /* CUDARNAINTERACTION_H_ */
