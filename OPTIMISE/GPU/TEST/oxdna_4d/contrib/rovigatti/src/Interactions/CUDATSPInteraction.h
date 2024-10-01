/*
 * CUDATSPInteraction.h
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#ifndef CUDATSPINTERACTION_H_
#define CUDATSPINTERACTION_H_

#define TSP_MAX_ARMS 20

#include "CUDA/Interactions/CUDABaseInteraction.h"
#include "TSPInteraction.h"

/**
 * @brief Contains information about anchors' bonded neighbours
 */
typedef struct
__align__(8) {
	int n[TSP_MAX_ARMS];
} TSP_anchor_bonds;

/**
 * @brief Handles interactions between TSPs on CUDA. See TSPInteraction for a list of options.
 */

class CUDATSPInteraction: public CUDABaseInteraction, public TSPInteraction {
protected:
	int *_h_anchors, *_d_anchors;
	TSP_anchor_bonds *_h_anchor_neighs, *_d_anchor_neighs;

	void _setup_anchors();
public:
	CUDATSPInteraction();
	virtual ~CUDATSPInteraction();

	c_number get_cuda_rcut() {
		return this->get_rcut();
	}
	void get_settings(input_file &inp);

	void cuda_init(int N);

	void compute_forces(CUDABaseList*lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox*d_box);
};

extern "C" BaseInteraction *make_CUDATSPInteraction() {
	return new CUDATSPInteraction();
}

#endif /* CUDATSPINTERACTION_H_ */
