/*
 * CUDABussiThermostat.cpp
 *
 *  Created on: Feb 15, 2013
 *      Author: rovigatti
 */

#include "CUDABussiThermostat.h"

#include <curand_kernel.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/transform_reduce.h>

template<typename number4>
struct sum_K {
    __device__ number4 operator()(const number4& a, const number4& b) const {
    	number4 res;
        res.x = res.y = res.z = 0;
        res.w = a.w + b.w;
        return res;
    }
 };

template<typename number4>
struct compute_K {
    __device__ number4 operator()(const number4& a) const {
    	number4 res;
    	res.x = res.y = res.z = 0;
        res.w = 0.5f * CUDA_DOT(a, a);
        return res;
    }
 };

template<typename number, typename number4>
__global__ void bussi_thermostat(number4 *vels, number4 *Ls, number rescale_factor_t, number rescale_factor_r, int N) {
	if(IND >= N) return;

	number4 v = vels[IND];
	v.x *= rescale_factor_t;
	v.y *= rescale_factor_t;
	v.z *= rescale_factor_t;
	v.w = (v.x*v.x + v.y*v.y + v.z*v.z) * (number) 0.5f;
	vels[IND] = v;

	number4 L = Ls[IND];
	L.x *= rescale_factor_r;
	L.y *= rescale_factor_r;
	L.z *= rescale_factor_r;
	L.w = (L.x*L.x + L.y*L.y + L.z*L.z) * (number) 0.5f;
	Ls[IND] = L;
}

template<typename number, typename number4>
CUDABussiThermostat<number, number4>::CUDABussiThermostat() : CUDABaseThermostat<number, number4>(), BussiThermostat<number>() {

}

template<typename number, typename number4>
CUDABussiThermostat<number, number4>::~CUDABussiThermostat() {

}

template<typename number, typename number4>
void CUDABussiThermostat<number, number4>::get_settings(input_file &inp) {
	BussiThermostat<number>::get_settings(inp);
	CUDABaseThermostat<number, number4>::get_cuda_settings(inp);
}

template<typename number, typename number4>
void CUDABussiThermostat<number, number4>::init(int N) {
	BussiThermostat<number>::init(N);

	this->_setup_rand(N);
}

template<typename number, typename number4>
bool CUDABussiThermostat<number, number4>::would_activate(llint curr_step) {
	return (curr_step % this->_newtonian_steps == 0);
}

template<typename number, typename number4>
void CUDABussiThermostat<number, number4>::apply_cuda(number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_vels, number4 *d_Ls, llint curr_step) {
	if(!would_activate(curr_step)) return;

	// we first calculate the current kinetic energy
	thrust::device_ptr<number4> t_vels = thrust::device_pointer_cast(d_vels);
	thrust::device_ptr<number4> t_Ls = thrust::device_pointer_cast(d_Ls);

	number4 zero = {0., 0., 0., 0.};
	number4 K_now_t = thrust::transform_reduce(t_vels, t_vels + this->_N_part, compute_K<number4>(), zero, sum_K<number4>());
	number4 K_now_r = thrust::transform_reduce(t_Ls, t_Ls + this->_N_part, compute_K<number4>(), zero, sum_K<number4>());

	this->_update_K(this->_K_t);
	this->_update_K(this->_K_r);

	number rescale_factor_t = sqrt(this->_K_t / K_now_t.w);
	number rescale_factor_r = sqrt(this->_K_r / K_now_r.w);
	//printf("%lf %lf %lf %lf\n", K_now_t.w, K_now_r.w, rescale_factor_t, rescale_factor_r);

	bussi_thermostat<number, number4>
		<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
		(d_vels, d_Ls, rescale_factor_t, rescale_factor_r, this->_N_part);
}

template class CUDABussiThermostat<float, float4>;
template class CUDABussiThermostat<double, LR_double4>;
