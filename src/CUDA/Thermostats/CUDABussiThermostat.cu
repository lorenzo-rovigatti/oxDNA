/*
 * CUDABussiThermostat.cpp
 *
 *  Created on: Feb 15, 2013
 *      Author: rovigatti
 */

#include "CUDABussiThermostat.h"

#include "../../Utilities/ConfigInfo.h"

#include <curand_kernel.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/transform_reduce.h>

struct sum_K {
	__device__
	c_number4 operator()(const c_number4& a, const c_number4& b) const {
		c_number4 res;
		res.x = res.y = res.z = 0;
		res.w = a.w + b.w;
		return res;
	}
};

struct compute_K {
	__device__
	c_number4 operator()(const c_number4& a) const {
		c_number4 res;
		res.x = res.y = res.z = 0;
		res.w = 0.5f * CUDA_DOT(a, a);
		return res;
	}
};

struct compute_K_without_com {
	c_number4 v_com;

	compute_K_without_com(c_number4 my_v_com) : v_com(my_v_com) {

	}

	__device__
	c_number4 operator()(const c_number4& v) const {
		c_number4 v_rel = v - v_com;
		c_number4 res;
		res.x = res.y = res.z = 0;
		res.w = 0.5f * CUDA_DOT(v_rel, v_rel);
		return res;
	}
};


__global__ void bussi_thermostat(c_number4 *vels, c_number4 *Ls, c_number4 v_com, c_number rescale_factor_t, c_number rescale_factor_r, int N) {
	if(IND >= N) return;

	c_number4 v = vels[IND];
	v.x = (v.x - v_com.x) * rescale_factor_t + v_com.x;
	v.y = (v.y - v_com.y) * rescale_factor_t + v_com.y;
	v.z = (v.z - v_com.z) * rescale_factor_t + v_com.z;
	v.w = (v.x * v.x + v.y * v.y + v.z * v.z) * (c_number) 0.5f;
	vels[IND] = v;

	c_number4 L = Ls[IND];
	L.x *= rescale_factor_r;
	L.y *= rescale_factor_r;
	L.z *= rescale_factor_r;
	L.w = (L.x * L.x + L.y * L.y + L.z * L.z) * (c_number) 0.5f;
	Ls[IND] = L;
}

CUDABussiThermostat::CUDABussiThermostat() :
				CUDABaseThermostat(),
				BussiThermostat() {

}

CUDABussiThermostat::~CUDABussiThermostat() {

}

void CUDABussiThermostat::get_settings(input_file &inp) {
	BussiThermostat::get_settings(inp);
	CUDABaseThermostat::get_cuda_settings(inp);
}

void CUDABussiThermostat::init() {
	BussiThermostat::init();

	_setup_rand(CONFIG_INFO->N());
}

bool CUDABussiThermostat::would_activate(llint curr_step) {
	return (curr_step % _newtonian_steps == 0);
}

void CUDABussiThermostat::apply_cuda(c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_vels, c_number4 *d_Ls, llint curr_step) {
	if(!would_activate(curr_step)) return;

	int N = CONFIG_INFO->N();

	// we first calculate the current kinetic energy
	thrust::device_ptr<c_number4> t_vels = thrust::device_pointer_cast(d_vels);
	thrust::device_ptr<c_number4> t_Ls = thrust::device_pointer_cast(d_Ls);

	c_number4 zero = { 0., 0., 0., 0. };
	c_number4 v_com = thrust::reduce(t_vels, t_vels + N, zero);
	v_com *= 1. /  N;
	c_number4 K_now_t = thrust::transform_reduce(t_vels, t_vels + N, compute_K_without_com(v_com), zero, sum_K());
	c_number4 K_now_r = thrust::transform_reduce(t_Ls, t_Ls + N, compute_K(), zero, sum_K());

	_update_K(_K_t, _current_translational_degrees_of_freedom());
	_update_K(_K_r, _current_rotational_degrees_of_freedom());

	c_number rescale_factor_t = std::sqrt(_K_t / K_now_t.w);
	c_number rescale_factor_r = std::sqrt(_K_r / K_now_r.w);

	bussi_thermostat
		<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
		(d_vels, d_Ls, v_com, rescale_factor_t, rescale_factor_r, N);
}
