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

class CUDATEPInteraction: public CUDABaseInteraction, public TEPInteraction {
protected:
	c_number *_d_kt_pref, *_d_kb1_pref, *_d_kb2_pref, *_d_xk_bending, *_d_xu_bending;
	c_number4 *_d_o_vects, *_d_w_vects;
	llint _steps;
public:
	CUDATEPInteraction();
	virtual ~CUDATEPInteraction();

	void get_settings(input_file &inp) override;
	void cuda_init(int N) override;
	c_number get_cuda_rcut() {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList*lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox*d_box);
};

#endif /* CUDATEPINTERACTION_H_ */
