/*
 * CUDATEPInteraction.h
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#ifndef CUDATEPINTERACTION_H_
#define CUDATEPINTERACTION_H_

#include "CUDABaseInteraction.h"
#include "../../Interactions/TEPInteraction.h"

/**
 * @brief CUDA implementation of the {@link TEPInteraction TEP interaction}.
 */
template<typename number, typename number4>
class CUDATEPInteraction: public CUDABaseInteraction<number, number4>, public TEPInteraction<number> {
protected:
	number *_d_kt_pref, *_d_kb1_pref, *_d_kb2_pref, *_d_xk_bending, *_d_xu_bending;
	number4 *_d_o_vects, *_d_w_vects;
	llint _steps;
public:
	CUDATEPInteraction();
	virtual ~CUDATEPInteraction();

	void get_settings(input_file &inp);
	void cuda_init(number box_side, int N);
	number get_cuda_rcut() { return this->get_rcut(); }

	void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box);
};

#endif /* CUDATEPINTERACTION_H_ */
