/*
 * CUDAStarrInteraction.h
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#ifndef CUDASTARRINTERACTION_H_
#define CUDASTARRINTERACTION_H_

#include "CUDA/Interactions/CUDABaseInteraction.h"
#include "StarrInteraction.h"

#define HUB_SIZE 4

typedef struct
	__align__(8) {
		int n[HUB_SIZE - 1];
	} hub_bonds;

	/**
	 * @brief Handles interactions between Starr tetramers on CUDA.
	 */

	class CUDAStarrInteraction: public CUDABaseInteraction, public StarrInteraction {
	protected:
		int _N_hubs;
		int *_d_hubs;
		int *_d_strand_ids;
		hub_bonds *_d_hub_neighs;
		c_number4 *_d_n3_forces, *_d_n5_forces;

		void _setup_strand_ids();
		void _setup_hubs();
	public:
		CUDAStarrInteraction();
		virtual ~CUDAStarrInteraction();

		c_number get_cuda_rcut() {
			return this->get_rcut();
		}
		void get_settings(input_file &inp);

		void cuda_init(int N);

		void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box);
	};

	extern "C" BaseInteraction *make_CUDAStarrInteraction();

#endif /* CUDASTARRINTERACTION_H_ */
