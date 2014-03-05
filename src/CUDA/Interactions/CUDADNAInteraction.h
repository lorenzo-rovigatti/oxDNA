/*
 * CUDADNAInteraction.h
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#ifndef CUDADNAINTERACTION_H_
#define CUDADNAINTERACTION_H_

#include "CUDABaseInteraction.h"
#include "../../Interactions/DNAInteraction.h"

/**
 * @brief CUDA implementation of the oxDNA model, as provided by DNAInteraction.
 */
template<typename number, typename number4>
class CUDADNAInteraction: public CUDABaseInteraction<number, number4>, public DNAInteraction<number> {
public:
	CUDADNAInteraction();
	virtual ~CUDADNAInteraction();

	void get_settings(input_file &inp);
	void cuda_init(number box_side, int N);
	number get_cuda_rcut() { return this->get_rcut(); }

	void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, LR_GPU_matrix<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds);
};

#endif /* CUDADNAINTERACTION_H_ */
