/*
 * CUDALangevinThermostat.cu
 *
 *  Created on: Sep 1, 2014
 *      Author: mzimmer
 */

#include <curand_kernel.h>
#include "CUDALangevinThermostat.h"

__global__ void langevin_thermostat(curandState *rand_state, c_number4 *vels, c_number4 *Ls, c_number _dt, c_number gamma_trans, c_number rescale_factor_trans, c_number gamma_rot, c_number rescale_factor_rot, int N) {
	if(IND < N) {
		curandState state = rand_state[IND];
		c_number4 vFuzz;
		c_number4 LFuzz;
		gaussian(state, vFuzz.x, vFuzz.y);
		gaussian(state, vFuzz.z, LFuzz.x);
		gaussian(state, LFuzz.y, LFuzz.z);

		c_number4 v = vels[IND];
		c_number4 L = Ls[IND];

		//Could define operators for GPU_quat
		v.x += _dt * (-gamma_trans * v.x + vFuzz.x * rescale_factor_trans);
		v.y += _dt * (-gamma_trans * v.y + vFuzz.y * rescale_factor_trans);
		v.z += _dt * (-gamma_trans * v.z + vFuzz.z * rescale_factor_trans);
		L.x += _dt * (-gamma_rot * L.x + LFuzz.x * rescale_factor_rot);
		L.y += _dt * (-gamma_rot * L.y + LFuzz.y * rescale_factor_rot);
		L.z += _dt * (-gamma_rot * L.z + LFuzz.z * rescale_factor_rot);

		//if (IND==5) printf("Diff is: %f %f %f \n",L.x-Ls[IND].x, L.y-Ls[IND].y, L.z-Ls[IND].z);

		vels[IND] = v;
		Ls[IND] = L;

		rand_state[IND] = state;
	}
}

CUDALangevinThermostat::CUDALangevinThermostat() :
				CUDABaseThermostat(),
				LangevinThermostat() {

}

CUDALangevinThermostat::~CUDALangevinThermostat() {

}

void CUDALangevinThermostat::get_settings(input_file &inp) {
	LangevinThermostat::get_settings(inp);
	CUDABaseThermostat::get_cuda_settings(inp);
}

void CUDALangevinThermostat::init(int N) {
	LangevinThermostat::init(N);

	this->_setup_rand(N);
}

bool CUDALangevinThermostat::would_activate(llint curr_step) {
	return 1;
}

void CUDALangevinThermostat::apply_cuda(c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_vels, c_number4 *d_Ls, llint curr_step) {
	if(!would_activate(curr_step)) return;

	langevin_thermostat
	<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
	(this->_d_rand_state, d_vels, d_Ls, this->_dt, this->_gamma_trans, this-> _rescale_factor_trans, this->_gamma_rot, this->_rescale_factor_rot, this->_N_part);
}
