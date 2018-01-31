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

typedef struct __align__(8) {
	int n[CENTRE_N_NEIGHS-1];
} centre_bonds;

/**
 * @brief Handles interactions between Levy tetramers on CUDA.
 */
template<typename number, typename number4>
class CUDALevyInteraction: public CUDABaseInteraction<number, number4>, public LevyInteraction<number> {
protected:
	int *_d_centres;
	centre_bonds *_d_centre_neighs;
	number4 *_d_n3_forces, *_d_n5_forces;

	void _setup_centres();
public:
	CUDALevyInteraction();
	virtual ~CUDALevyInteraction();

	number get_cuda_rcut() { return this->get_rcut(); }
	void get_settings(input_file &inp);

	void cuda_init(number box_side, int N);

	void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box);
};

extern "C" IBaseInteraction<float> *make_CUDALevyInteraction_float();
extern "C" IBaseInteraction<double> *make_CUDALevyInteraction_double();

#endif /* CUDALEVYINTERACTION_H_ */
