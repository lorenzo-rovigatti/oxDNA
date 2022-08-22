/*
 * CUDABaseInteraction.h
 *
 *  Created on: 18/feb/2013
 *      Author: lorenzo
 */

#ifndef CUDABASEINTERACTION_H_
#define CUDABASEINTERACTION_H_

#include "../CUDAUtils.h"
#include "../Lists/CUDABaseList.h"
#include "../cuda_utils/CUDABox.h"

/**
 * @brief Abstract class providing an interface for CUDA-based interactions.
 */

class CUDABaseInteraction {
protected:
	bool _use_edge;
	/// c_number of force slots per particle. Used only if _use_edge == true.
	int _n_forces;
	CUDA_kernel_cfg _launch_cfg;
	CUDA_kernel_cfg _ffs_hb_precalc_kernel_cfg;
	CUDA_kernel_cfg _ffs_dist_precalc_kernel_cfg;
	CUDA_kernel_cfg _ffs_hb_eval_kernel_cfg;
	CUDA_kernel_cfg _ffs_dist_eval_kernel_cfg;

	int _N;
	c_number _box_side;

	c_number4 *_d_edge_forces;
	c_number4 *_d_edge_torques;

	virtual void _sum_edge_forces(c_number4 *d_forces);
	virtual void _sum_edge_forces_torques(c_number4 *d_forces, c_number4 *d_torques);
public:
	CUDABaseInteraction();
	virtual ~CUDABaseInteraction();

	virtual void get_settings(input_file &inp) = 0;
	virtual void get_cuda_settings(input_file &inp);
	virtual void cuda_init(c_number box_side, int N);
	virtual c_number get_cuda_rcut() = 0;

	virtual void sync_host() {}
	virtual void sync_GPU() {}

	void set_launch_cfg(CUDA_kernel_cfg &launch_cfg);

	virtual void compute_forces(CUDABaseList*lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox*d_box) = 0;
	virtual void _hb_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, float *hb_energies, int n_threads, bool *region_is_nearhb, CUDA_kernel_cfg hb_kernel_cfg, CUDABox*d_box);
	virtual void _near_hb_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, bool *nearly_bonded_array, int n_threads, bool *region_is_nearhb, CUDA_kernel_cfg hb_kernel_cfg, CUDABox*d_box);
	virtual void _dist_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, c_number *op_dists, int n_threads, CUDA_kernel_cfg dist_kernel_cfg, CUDABox*d_box);
};

#endif /* CUDABASEINTERACTION_H_ */
