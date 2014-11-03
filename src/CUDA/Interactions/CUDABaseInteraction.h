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

/**
 * @brief Abstract class providing an interface for CUDA-based interactions.
 */
template<typename number, typename number4>
class CUDABaseInteraction {
protected:
	bool _use_edge;
	/// Number of force slots per particle. Used only if _use_edge == true.
	int _n_forces;
	CUDA_kernel_cfg _launch_cfg;
	CUDA_kernel_cfg _ffs_hb_precalc_kernel_cfg;
	CUDA_kernel_cfg _ffs_dist_precalc_kernel_cfg;
	CUDA_kernel_cfg _ffs_hb_eval_kernel_cfg;
	CUDA_kernel_cfg _ffs_dist_eval_kernel_cfg;

	int _N;
	number _box_side;

	number4 *_d_edge_forces;
	number4 *_d_edge_torques;

	virtual void _sum_edge_forces(number4 *d_forces);
	virtual void _sum_edge_forces_torques(number4 *d_forces, number4 *d_torques);
public:
	CUDABaseInteraction();
	virtual ~CUDABaseInteraction();

	virtual void get_settings(input_file &inp) = 0;
	virtual void get_cuda_settings(input_file &inp);
	virtual void cuda_init(number box_side, int N);
	virtual number get_cuda_rcut() = 0;

	void set_launch_cfg(CUDA_kernel_cfg &launch_cfg);

	virtual void compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_qorientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds) = 0;
	virtual void _hb_op_precalc(number4 *poss, GPU_quat<number> *orientations, int *op_pairs1, int *op_pairs2, float *hb_energies, int n_threads, bool *region_is_nearhb, CUDA_kernel_cfg hb_kernel_cfg);
	virtual void _near_hb_op_precalc(number4 *poss, GPU_quat<number> *orientations, int *op_pairs1, int *op_pairs2, bool *nearly_bonded_array, int n_threads, bool *region_is_nearhb, CUDA_kernel_cfg hb_kernel_cfg);
	virtual void _dist_op_precalc(number4 *poss, GPU_quat<number> *orientations, int *op_pairs1, int *op_pairs2, number *op_dists, int n_threads, CUDA_kernel_cfg dist_kernel_cfg);
};

#endif /* CUDABASEINTERACTION_H_ */
