/*
 * CUDAStarrInteraction.h
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#ifndef CUDASTARRINTERACTION_H_
#define CUDASTARRINTERACTION_H_

#include "CUDABaseInteraction.h"
#include "StarrInteraction.h"

#define HUB_SIZE 4

typedef struct __align__(8) {
	int n[HUB_SIZE-1];
} hub_bonds;

/**
 * @brief Handles interactions between Starr tetramers on CUDA.
 */
template<typename number, typename number4>
class CUDAStarrInteraction: public CUDABaseInteraction<number, number4>, public StarrInteraction<number> {
protected:
	int _N_hubs;
	int *_d_hubs;
	int *_d_strand_ids;
	hub_bonds *_d_hub_neighs;
	number4 *_d_n3_forces, *_d_n5_forces;

	void _setup_strand_ids();
	void _setup_hubs();
public:
	CUDAStarrInteraction();
	virtual ~CUDAStarrInteraction();

	number get_cuda_rcut() { return this->get_rcut(); }
	void get_settings(input_file &inp);

	void cuda_init(number box_side, int N);

	void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds);
};

extern "C" IBaseInteraction<float> *make_CUDAStarrInteraction_float();
extern "C" IBaseInteraction<double> *make_CUDAStarrInteraction_double();

#endif /* CUDASTARRINTERACTION_H_ */
