/**
 * @file    CUDABaseBackend.h
 * @date    25/nov/2010
 * @author  lorenzo
 *
 *
 */

#ifndef CUDABASEBACKEND_H_
#define CUDABASEBACKEND_H_

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "../cuda_utils/cuda_device_utils.h"
#include "../CUDAUtils.h"
#include "../Lists/CUDABaseList.h"
#include "../CUDA_sort.cuh"
#include "../Interactions/CUDABaseInteraction.h"
#include "../../Observables/BaseObservable.h"

/**
 * @brief Basic simulation backend on CUDA. All CUDA backends should inherit from this class as well as from a regular CPU backend
 *
 * This class does not actually do any computation but provides basic CUDA facilities.
 *
 * @verbatim
[CUDA_device = <int> (CUDA-enabled device to run the simulation on. If it is not specified or it is given a negative number, a suitable device will be automatically chosen.)]
[CUDA_sort_every = <int> (sort particles according to a 3D Hilbert curve every CUDA_sort_every time steps. This will greatly enhnance performances for some types of interaction. Defaults to 0, which disables sorting.)]
[threads_per_block = <int> (Number of threads per block on the CUDA grid. defaults to 2 * the size of a warp.)]
@endverbatim
 */
template<typename number, typename number4>
class CUDABaseBackend {
private:
	ConfigInfo<number> _prv_config_info;

protected:
	/// if 0 then do not sort. If it's > 1 then sort particles every _sort_every updates
	int _sort_every;
	int _device_number;
	cudaDeviceProp _device_prop;
	CUDA_kernel_cfg _particles_kernel_cfg;
	CUDABaseList<number, number4> *_cuda_lists;
	size_t _vec_size, _bonds_size, _orient_size;
	number _sqr_verlet_skin;

	/// used for sorting
	number4 *_d_buff_poss;
	GPU_quat<number> *_d_buff_orientations;
	LR_bonds *_d_buff_bonds;
	int *_d_hindex, *_d_sorted_hindex, *_d_inv_sorted_hindex;

	number4 *_d_poss, *_h_poss;
	LR_bonds *_d_bonds, *_h_bonds;
	GPU_quat<number> *_d_orientations, *_h_orientations;
	number4 *_d_list_poss;
	/// It is stored in pinned memory, i.e. on the host but it can be accessed directly from the device
	bool *_d_are_lists_old;

	CUDABaseInteraction<number, number4> *_cuda_interaction;

	virtual void _host_to_gpu();
	virtual void _gpu_to_host();

	virtual void _sort_index();
	virtual void _init_CUDA_kernel_cfgs();

	/**
	 * @brief This handy method automatically selects a CUDA device to run the simulation on if there are no user-specified devices.
	 */
	virtual void _choose_device();
public:
	CUDABaseBackend();
	virtual ~CUDABaseBackend();

	virtual void get_settings(input_file &inp);
	virtual void init_cuda(ConfigInfo<number> &config_info);
};

#endif /* CUDABASEBACKEND_H_ */
