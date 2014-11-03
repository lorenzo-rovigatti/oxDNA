/*
 * CUDATSPInteraction.h
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#ifndef CUDATSPINTERACTION_H_
#define CUDATSPINTERACTION_H_

#define TSP_MAX_ARMS 20

#include "CUDABaseInteraction.h"
#include "../../Interactions/TSPInteraction.h"

/**
 * @brief Contains information about anchors' bonded neighbours
 */
typedef struct __align__(8) {
	int n[TSP_MAX_ARMS];
} TSP_anchor_bonds;

/**
 * @brief Handles interactions between TSPs on CUDA. See TSPInteraction for a list of options.
 */
template<typename number, typename number4>
class CUDATSPInteraction: public CUDABaseInteraction<number, number4>, public TSPInteraction<number> {
protected:
	int *_h_anchors, *_d_anchors;
	TSP_anchor_bonds *_h_anchor_neighs, *_d_anchor_neighs;

	void _setup_anchors();
public:
	CUDATSPInteraction();
	virtual ~CUDATSPInteraction();

	number get_cuda_rcut() { return this->get_rcut(); }
	void get_settings(input_file &inp);

	void cuda_init(number box_side, int N);

	void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds);
};

#endif /* CUDATSPINTERACTION_H_ */
