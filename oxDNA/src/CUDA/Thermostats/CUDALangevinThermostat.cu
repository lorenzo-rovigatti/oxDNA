/*
 * CUDALangevinThermostat.cu
 *
 *  Created on: Sep 1, 2014
 *      Author: mzimmer
 */

#include <curand_kernel.h>
#include "CUDALangevinThermostat.h"

template<typename number, typename number4>
__global__ void langevin_thermostat(curandState *rand_state, number4 *vels, number4 *Ls, number _dt, number gamma_trans, number rescale_factor_trans, number gamma_rot, number rescale_factor_rot, int N) {
	if (IND < N) {
		curandState state = rand_state[IND];
		number4 vFuzz;
		number4 LFuzz;		
		gaussian(state, vFuzz.x, vFuzz.y);
		gaussian(state, vFuzz.z, LFuzz.x);
		gaussian(state, LFuzz.y, LFuzz.z);
		
		number4 v = vels[IND];
		number4 L = Ls[IND];

		//Could define operators for GPU_quat
		v.x += _dt * (-gamma_trans*v.x + vFuzz.x*rescale_factor_trans);
		v.y += _dt * (-gamma_trans*v.y + vFuzz.y*rescale_factor_trans);
		v.z += _dt * (-gamma_trans*v.z + vFuzz.z*rescale_factor_trans);
		L.x += _dt * (-gamma_rot*L.x + LFuzz.x*rescale_factor_rot);
		L.y += _dt * (-gamma_rot*L.y + LFuzz.y*rescale_factor_rot);
		L.z += _dt * (-gamma_rot*L.z + LFuzz.z*rescale_factor_rot);

		//if (IND==5) printf("Diff is: %f %f %f \n",L.x-Ls[IND].x, L.y-Ls[IND].y, L.z-Ls[IND].z);

		vels[IND] = v;
		Ls[IND] = L;

		rand_state[IND] = state;
	}
}

template<typename number, typename number4>
CUDALangevinThermostat<number, number4>::CUDALangevinThermostat() : CUDABaseThermostat<number, number4>(), LangevinThermostat<number>() {

}

template<typename number, typename number4>
CUDALangevinThermostat<number, number4>::~CUDALangevinThermostat() {

}

template<typename number, typename number4>
void CUDALangevinThermostat<number, number4>::get_settings(input_file &inp) {
	LangevinThermostat<number>::get_settings(inp);
	CUDABaseThermostat<number, number4>::get_cuda_settings(inp);
}

template<typename number, typename number4>
void CUDALangevinThermostat<number, number4>::init(int N) {
	LangevinThermostat<number>::init(N);

	this->_setup_rand(N);
}

template<typename number, typename number4>
bool CUDALangevinThermostat<number, number4>::would_activate(llint curr_step) {
	return 1;
}

template<typename number, typename number4>
void CUDALangevinThermostat<number, number4>::apply_cuda(number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_vels, number4 *d_Ls, llint curr_step) {
	if(!would_activate(curr_step)) return;
	//printf("%d \n",this->_N_part);
	langevin_thermostat<number, number4>
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(this->_d_rand_state, d_vels, d_Ls, this->_dt, this->_gamma_trans, this-> _rescale_factor_trans, this->_gamma_rot, this->_rescale_factor_rot, this->_N_part);
}

template class CUDALangevinThermostat<float, float4>;
template class CUDALangevinThermostat<double, LR_double4>;
