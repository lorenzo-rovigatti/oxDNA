/*
 * CUDALevyInteraction.h
 *
 *  Created on: 29/jun/2016
 *      Author: lorenzo
 */

#ifndef CUDALEVYINTERACTION_H_
#define CUDALEVYINTERACTION_H_

#include "CUDA/Interactions/CUDABaseInteraction.h"
#include "LevyInteraction.h"

#define CENTRE_N_NEIGHS 4

typedef struct
	__align__(8) {
		int n[CENTRE_N_NEIGHS - 1];
	} centre_bonds;

	/**
	 * @brief Handles interactions between Levy tetramers on CUDA.
	 */

	class CUDALevyInteraction: public CUDABaseInteraction, public LevyInteraction {
	protected:
		int *_d_centres;
		centre_bonds *_d_centre_neighs;
		c_number4 *_d_n3_forces, *_d_n5_forces;

		void _setup_centres();
	public:
		CUDALevyInteraction();
		virtual ~CUDALevyInteraction();

		c_number get_cuda_rcut() {
			return this->get_rcut();
		}
		void get_settings(input_file &inp);

		void cuda_init(int N);

		void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box);
	};

	extern "C" BaseInteraction *make_CUDALevyInteraction();

#endif /* CUDALEVYINTERACTION_H_ */
