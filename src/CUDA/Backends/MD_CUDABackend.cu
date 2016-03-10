/*
 * MD_CUDABackend.cpp
 *
 *  Created on: 03/set/2010
 *      Author: lorenzo
 */

#include "MD_CUDABackend.h"
#include "CUDA_MD.cuh"
#include "../CUDA_base_interactions.h"
#include "../../Interactions/DNAInteraction.h"
#include "../Thermostats/CUDAThermostatFactory.h"
#include <thrust/sort.h>
#include <typeinfo>

#include "../../Forces/COMForce.h"
#include "../../Forces/ConstantRateForce.h"
#include "../../Forces/ConstantRateTorque.h"
#include "../../Forces/ConstantTrap.h"
#include "../../Forces/LowdimMovingTrap.h"
#include "../../Forces/MovingTrap.h"
#include "../../Forces/MutualTrap.h"
#include "../../Forces/RepulsionPlane.h"
#include "../../Forces/RepulsionPlaneMoving.h"
#include "../../Forces/RepulsiveSphere.h"
#include "../../Forces/LJWall.h"

template<typename number, typename number4>
MD_CUDABackend<number, number4>::MD_CUDABackend() : MDBackend<number>(), CUDABaseBackend<number, number4>(), _max_ext_forces(0) {
	this->_is_CUDA_sim = true;
	_use_edge = false;
	_any_rigid_body = false;

	_d_vels = _d_Ls = _d_forces = _d_torques = NULL;
	_h_vels = _h_Ls = _h_forces = _h_torques = NULL;
	_h_gpu_index = _h_cpu_index = NULL;

	_cuda_thermostat = NULL;

	_h_ext_forces = NULL;
	_d_ext_forces = NULL;

	_restart_step_counter = false;
}

template<typename number, typename number4>
MD_CUDABackend<number, number4>::~MD_CUDABackend() {
	if(_d_vels != NULL) {
		CUDA_SAFE_CALL( cudaFree(_d_vels) );
		CUDA_SAFE_CALL( cudaFree(_d_Ls) );
		CUDA_SAFE_CALL( cudaFree(_d_forces) );
		CUDA_SAFE_CALL( cudaFree(_d_torques) );
	}

	if(this->_sort_every > 0) {
		CUDA_SAFE_CALL( cudaFree(_d_buff_vels) );
		CUDA_SAFE_CALL( cudaFree(_d_buff_Ls) );
	}

	if(_h_gpu_index != NULL) {
		delete[] _h_gpu_index;
		delete[] _h_cpu_index;
	}

	if(this->_external_forces) {
		if(_h_ext_forces != NULL)
			delete[] _h_ext_forces;
		if(_d_ext_forces != NULL)
			CUDA_SAFE_CALL( cudaFree(_d_ext_forces) );
	}

	if(_h_vels != NULL) {
		delete[] _h_vels;
		delete[] _h_Ls;
		delete[] _h_forces;
		delete[] _h_torques;
	}

	if(_cuda_thermostat != NULL) delete _cuda_thermostat;
}

template<typename number, typename number4>
void MD_CUDABackend<number, number4>::_host_to_gpu() {
	CUDABaseBackend<number, number4>::_host_to_gpu();
	CUDA_SAFE_CALL( cudaMemcpy(_d_vels, _h_vels, this->_vec_size, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_Ls, _h_Ls, this->_vec_size, cudaMemcpyHostToDevice) );
}

template<typename number, typename number4>
void MD_CUDABackend<number, number4>::_gpu_to_host() {
	CUDABaseBackend<number, number4>::_gpu_to_host();
	CUDA_SAFE_CALL( cudaMemcpy(_h_vels, _d_vels, this->_vec_size, cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(_h_Ls, _d_Ls, this->_vec_size, cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(_h_forces, _d_forces, this->_vec_size, cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(_h_torques, _d_torques, this->_vec_size, cudaMemcpyDeviceToHost) );
}

template<typename number, typename number4>
void MD_CUDABackend<number, number4>::_host_particles_to_gpu() {
	for(int i = 0; i < this->_N; i++) {
		int gpu_index = _h_gpu_index[i];
		BaseParticle<number> *p = this->_particles[gpu_index];

		this->_h_poss[i].x = p->pos.x;
		this->_h_poss[i].y = p->pos.y;
		this->_h_poss[i].z = p->pos.z;

		// convert index and type into a float
		int msk = -1 << 22;// bynary mask: all 1's and 22 0's;
		// btype has a sign, and thus has to go first
		this->_h_poss[i].w = GpuUtils::int_as_float( (p->btype << 22) | ((~msk) & p->index) );
		// we immediately check that the index and base type that we read are sensible
		int mybtype = (GpuUtils::float_as_int(this->_h_poss[i].w)) >> 22;
		int myindex = (GpuUtils::float_as_int(this->_h_poss[i].w)) & (~msk);
		if (p->btype != mybtype) {
			throw oxDNAException ("Could not treat the type (A, C, G, T or something specific) of particle %d; On CUDA, the maximum \"unique\" identity is 512");
		}
		if (p->index != myindex) {
			throw oxDNAException ("Could not treat the index of particle %d; remember that on CUDA the maximum number of particles is 2^21", p->index);
		}

		if (p->n3 == P_VIRTUAL) this->_h_bonds[i].n3 = P_INVALID;
		else this->_h_bonds[i].n3 = _h_cpu_index[p->n3->index];
		if (p->n5 == P_VIRTUAL) this->_h_bonds[i].n5 = P_INVALID;
		else this->_h_bonds[i].n5 = _h_cpu_index[p->n5->index];

		_h_vels[i].x = p->vel.x;
		_h_vels[i].y = p->vel.y;
		_h_vels[i].z = p->vel.z;

		_h_Ls[i].x = p->L.x;
		_h_Ls[i].y = p->L.y;
		_h_Ls[i].z = p->L.z;
	
		number trace = p->orientation.v1.x+p->orientation.v2.y+p->orientation.v3.z;
		if(trace > 0) {
			number s = .5/sqrt(trace + 1);
			this->_h_orientations[i].w = .25/s;
			this->_h_orientations[i].x = (p->orientation.v3.y-p->orientation.v2.z)*s;
			this->_h_orientations[i].y = (p->orientation.v1.z-p->orientation.v3.x)*s;
			this->_h_orientations[i].z = (p->orientation.v2.x-p->orientation.v1.y)*s;
		} else {    //Finding largest diagonal element
			if ( (p->orientation.v1.x > p->orientation.v2.y) && (p->orientation.v1.x > p->orientation.v3.z) ) { 
				number s = .5/sqrt(1+p->orientation.v1.x-p->orientation.v2.y-p->orientation.v3.z);
				this->_h_orientations[i].w = (p->orientation.v3.y-p->orientation.v2.z)*s;
				this->_h_orientations[i].x = .25/s;
				this->_h_orientations[i].y = (p->orientation.v1.y+p->orientation.v2.x)*s;
				this->_h_orientations[i].z = (p->orientation.v1.z+p->orientation.v3.x)*s;
			} else if (p->orientation.v2.y > p->orientation.v3.z) {
				number s = .5/sqrt(1+p->orientation.v2.y-p->orientation.v1.x-p->orientation.v3.z);
				this->_h_orientations[i].w = (p->orientation.v1.z-p->orientation.v3.x)*s;
				this->_h_orientations[i].x = (p->orientation.v1.y+p->orientation.v2.x)*s;
				this->_h_orientations[i].y = .25/s;
				this->_h_orientations[i].z = (p->orientation.v2.z+p->orientation.v3.y)*s;
			} else {
				number s = .5/sqrt(1+p->orientation.v3.z-p->orientation.v1.x-p->orientation.v2.y);
				this->_h_orientations[i].w = (p->orientation.v2.x-p->orientation.v1.y)*s;
				this->_h_orientations[i].x = (p->orientation.v1.z+p->orientation.v3.x)*s;
				this->_h_orientations[i].y = (p->orientation.v2.z+p->orientation.v3.y)*s;
				this->_h_orientations[i].z = .25/s;
			}
		}
	} 
	this->_host_to_gpu();
}

template<typename number, typename number4>
void MD_CUDABackend<number, number4>::_gpu_to_host_particles() {
	this->_gpu_to_host();

	for(int i = 0; i < this->_N; i++) {
		// since we may have been sorted all the particles in a different order
		// we first take the particle index from the 4th component of its
		// position, and then use that index to access the right BaseParticle pointer
		int msk = (-1 << 22);
		int newindex = ((GpuUtils::float_as_int(this->_h_poss[i].w)) & (~msk));
		_h_gpu_index[i] = newindex;
		_h_cpu_index[newindex] = i;
		BaseParticle<number> *p = this->_particles[newindex];
		assert(p->index == newindex);

		p->pos.x = this->_h_poss[i].x;
		p->pos.y = this->_h_poss[i].y;
		p->pos.z = this->_h_poss[i].z;
		// get index and type for the fourth component of the position
		p->btype = (GpuUtils::float_as_int(this->_h_poss[i].w)) >> 22;

		if (this->_h_bonds[i].n3 == P_INVALID) p->n3 = P_VIRTUAL;
		else {
			int n3index = ((GpuUtils::float_as_int(this->_h_poss[this->_h_bonds[i].n3].w)) & (~msk));
			p->n3 = this->_particles[n3index];
		}
		if (this->_h_bonds[i].n5 == P_INVALID) p->n5 = P_VIRTUAL;
		else {
			int n5index = ((GpuUtils::float_as_int(this->_h_poss[this->_h_bonds[i].n5].w)) & (~msk));
			p->n5 = this->_particles[n5index];
		}

		p->vel.x = _h_vels[i].x;
		p->vel.y = _h_vels[i].y;
		p->vel.z = _h_vels[i].z;

		p->L.x = _h_Ls[i].x;
		p->L.y = _h_Ls[i].y;
		p->L.z = _h_Ls[i].z;

		number sqx = this->_h_orientations[i].x*this->_h_orientations[i].x;
		number sqy = this->_h_orientations[i].y*this->_h_orientations[i].y;
		number sqz = this->_h_orientations[i].z*this->_h_orientations[i].z;
		number sqw = this->_h_orientations[i].w*this->_h_orientations[i].w;
		number xy = this->_h_orientations[i].x*this->_h_orientations[i].y;
		number xz = this->_h_orientations[i].x*this->_h_orientations[i].z;
		number xw = this->_h_orientations[i].x*this->_h_orientations[i].w;
		number yz = this->_h_orientations[i].y*this->_h_orientations[i].z;
		number yw = this->_h_orientations[i].y*this->_h_orientations[i].w;
		number zw = this->_h_orientations[i].z*this->_h_orientations[i].w;	
		number invs = 1 / (sqx + sqy + sqz + sqw);

		p->orientation.v1.x = (sqx-sqy-sqz+sqw)*invs;
		p->orientation.v1.y = 2*(xy-zw)*invs;
		p->orientation.v1.z = 2*(xz+yw)*invs;
		p->orientation.v2.x = 2*(xy+zw)*invs;
		p->orientation.v2.y = (-sqx+sqy-sqz+sqw)*invs;
		p->orientation.v2.z = 2*(yz-xw)*invs;
		p->orientation.v3.x = 2*(xz-yw)*invs;
		p->orientation.v3.y = 2*(yz+xw)*invs;
		p->orientation.v3.z = (-sqx-sqy+sqz+sqw)*invs;

		p->set_positions();
		p->orientationT = p->orientation.get_transpose();
	}
}

template<typename number, typename number4>
void MD_CUDABackend<number, number4>::_init_CUDA_MD_symbols() {
	float f_copy = this->_sqr_verlet_skin;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_sqr_verlet_skin, &f_copy, sizeof(float)) );
	f_copy = this->_dt;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_dt, &f_copy, sizeof(float)) );
	f_copy = this->_box_side;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_box_side, &f_copy, sizeof(float)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &this->_N, sizeof(int)) );
}

template<typename number, typename number4>
void MD_CUDABackend<number, number4>::_first_step() {
	first_step<number, number4>
		<<<this->_particles_kernel_cfg.blocks, this->_particles_kernel_cfg.threads_per_block>>>
		(this->_d_poss, this->_d_orientations, this->_d_list_poss, _d_vels, _d_Ls, _d_forces, _d_torques, this->_d_are_lists_old);
	CUT_CHECK_ERROR("_first_step error");
}

template<typename number, typename number4>
void MD_CUDABackend<number, number4>::_forces_second_step() {
	_set_external_forces();
	this->_cuda_interaction->compute_forces(this->_cuda_lists, this->_d_poss, this->_d_orientations, _d_forces, _d_torques, this->_d_bonds);

	second_step<number, number4>
		<<<this->_particles_kernel_cfg.blocks, this->_particles_kernel_cfg.threads_per_block>>>
		(this->_d_vels, this->_d_Ls, this->_d_forces, this->_d_torques);
		CUT_CHECK_ERROR("second_step");
}

template<typename number, typename number4>
void MD_CUDABackend<number, number4>::_set_external_forces() {
	set_external_forces<number, number4>
		<<<this->_particles_kernel_cfg.blocks, this->_particles_kernel_cfg.threads_per_block>>>
		(this->_d_poss, this->_d_orientations, _d_ext_forces, _d_forces, _d_torques, _curr_step, _max_ext_forces);
	CUT_CHECK_ERROR("set_external_forces");
}

template<typename number, typename number4>
void MD_CUDABackend<number, number4>::_sort_particles() {
	CUDABaseBackend<number, number4>::_sort_index();
	permute_particles<number, number4>
		<<<this->_particles_kernel_cfg.blocks, this->_particles_kernel_cfg.threads_per_block>>>
		(this->_d_sorted_hindex, this->_d_inv_sorted_hindex, this->_d_poss, _d_vels, _d_Ls, this->_d_orientations, this->_d_bonds, this->_d_buff_poss, _d_buff_vels, _d_buff_Ls, 				this->_d_buff_orientations, this->_d_buff_bonds);
		CUT_CHECK_ERROR("_permute_particles error");
		CUDA_SAFE_CALL( cudaMemcpy(this->_d_orientations, this->_d_buff_orientations, this->_orient_size, cudaMemcpyDeviceToDevice) );

	// copy back the sorted vectors
	CUDA_SAFE_CALL( cudaMemcpy(this->_d_poss, this->_d_buff_poss, this->_vec_size, cudaMemcpyDeviceToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(this->_d_bonds, this->_d_buff_bonds, this->_bonds_size, cudaMemcpyDeviceToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_vels, _d_buff_vels, this->_vec_size, cudaMemcpyDeviceToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_Ls, _d_buff_Ls, this->_vec_size, cudaMemcpyDeviceToDevice) );
}

template<typename number, typename number4>
void MD_CUDABackend<number, number4>::_thermalize(llint curr_step) {
	_cuda_thermostat->apply_cuda(this->_d_poss, this->_d_orientations, _d_vels, _d_Ls, curr_step);
}

template<typename number, typename number4>
void MD_CUDABackend<number, number4>::sim_step(llint curr_step) {
	this->_mytimer->resume();
	_curr_step = curr_step;

	this->_timer_first_step->resume();
	_first_step();
	cudaThreadSynchronize();
	this->_timer_first_step->pause();

	_timer_sorting->resume();
	if(this->_d_are_lists_old[0] && this->_sort_every > 0 && (this->_N_updates % this->_sort_every == 0)) {
		_sort_particles();
		cudaThreadSynchronize();
	}
	_timer_sorting->pause();

	this->_timer_lists->resume();
	if(this->_d_are_lists_old[0]) {
		this->_cuda_lists->update(this->_d_poss, this->_d_list_poss, this->_d_bonds);
		this->_d_are_lists_old[0] = false;
		this->_N_updates++;
		cudaThreadSynchronize();
	}
	this->_timer_lists->pause();

	this->_timer_forces->resume();
	_forces_second_step();
	cudaThreadSynchronize();
	this->_timer_forces->pause();

	this->_timer_thermostat->resume();
	_thermalize(curr_step);
	cudaThreadSynchronize();
	this->_timer_thermostat->pause();

	this->_mytimer->pause();
}

template<typename number, typename number4>
void MD_CUDABackend<number, number4>::get_settings(input_file &inp) {
	MDBackend<number>::get_settings(inp);
	CUDABaseBackend<number, number4>::get_settings(inp);
	int tmpi;
	if(getInputBoolAsInt(&inp, "use_edge", &tmpi, 0) == KEY_FOUND) {
		_use_edge = (tmpi != 0);
		if(_use_edge && sizeof(number) == sizeof(double)) throw oxDNAException("use_edge and double precision are not compatible");
	}
	
	getInputBool(&inp, "restart_step_counter", &_restart_step_counter, 1);

	_cuda_thermostat = CUDAThermostatFactory::make_thermostat<number, number4>(inp, this->_box_side);
	_cuda_thermostat->get_settings(inp);
}

template<typename number, typename number4>
void MD_CUDABackend<number, number4>::init(){
	MDBackend<number>::init();
	CUDABaseBackend<number, number4>::init_cuda(this->_config_info);

	_timer_sorting = TimingManager::instance()->new_timer(std::string("Hilbert sorting"), std::string("SimBackend"));

	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number4>(&_d_vels, this->_vec_size) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number4>(&_d_Ls, this->_vec_size) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number4>(&_d_forces, this->_vec_size) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number4>(&_d_torques, this->_vec_size) );

	CUDA_SAFE_CALL( cudaMemset(_d_forces, 0, this->_vec_size) );
	CUDA_SAFE_CALL( cudaMemset(_d_torques, 0, this->_vec_size) );

	_h_vels = new number4[this->_N];
	_h_Ls = new number4[this->_N];
	_h_forces = new number4[this->_N];
	_h_torques = new number4[this->_N];

	if(this->_external_forces) {
		if(this->_sort_every > 0) throw oxDNAException("External forces and CUDA_sort_every > 0 are not compatible");
		_h_ext_forces = new CUDA_trap<number>[this->_N * MAX_EXT_FORCES];
		CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<CUDA_trap<number> >(&_d_ext_forces, this->_N * MAX_EXT_FORCES * sizeof(CUDA_trap<number>)) );

		for(int i = 0; i < this->_N*MAX_EXT_FORCES; i++) _h_ext_forces[i].type = -1;

		ConstantRateForce<number> const_force;
		MutualTrap<number> mutual_trap;
		MovingTrap<number> moving_trap;
		LowdimMovingTrap<number> lowdim_moving_trap;
		RepulsionPlane<number> repulsion_plane;
		RepulsionPlaneMoving<number> repulsion_plane_moving;
		LJWall<number> LJ_wall;

		for(int i = 0; i < this->_N; i++) {
			BaseParticle<number> *p = this->_particles[i];

			for(int j = 0; j < p->N_ext_forces; j++) {
				_max_ext_forces = max(_max_ext_forces, p->N_ext_forces);

				CUDA_trap<number> *force = &(_h_ext_forces[j*this->_N + i]);
				if (typeid (*(p->ext_forces[j])) == typeid (const_force) ) {
					ConstantRateForce<number> * p_force = (ConstantRateForce<number> *) p->ext_forces[j];
					force->type = CUDA_TRAP_CONSTANT;
					force->constant.F0 = p_force->_F0;
					force->constant.rate = p_force->_rate;
					force->constant.x = p_force->_direction.x;
					force->constant.y = p_force->_direction.y;
					force->constant.z = p_force->_direction.z;
				}
				else if (typeid (*(p->ext_forces[j])) == typeid (mutual_trap) ) {
					MutualTrap<number> * p_force = (MutualTrap<number> *) p->ext_forces[j];
					force->type = CUDA_TRAP_MUTUAL;
					force->mutual.stiff = p_force->_stiff;
					force->mutual.r0 = p_force->_r0;
					force->mutual.p_ind = p_force->_p_ptr->index;
					force->mutual.PBC = p_force->PBC;
				}
				else if (typeid (*(p->ext_forces[j])) == typeid (moving_trap) ) {
					MovingTrap<number> * p_force = (MovingTrap<number> *) p->ext_forces[j];
					force->type = CUDA_TRAP_MOVING;
					force->moving.stiff = p_force->_stiff;
					force->moving.rate = p_force->_rate;
					force->moving.pos0 = make_float3 (p_force->_pos0.x, p_force->_pos0.y, p_force->_pos0.z);
					force->moving.dir = make_float3 (p_force->_direction.x, p_force->_direction.y, p_force->_direction.z);
				}
				else if (typeid (*(p->ext_forces[j])) == typeid (lowdim_moving_trap) ) {
					LowdimMovingTrap<number> * p_force = (LowdimMovingTrap<number> *) p->ext_forces[j];
					force->type = CUDA_TRAP_MOVING_LOWDIM;
					force->lowdim.stiff = p_force->_stiff;
					force->lowdim.rate = p_force->_rate;
					force->lowdim.pos0 = make_float3 (p_force->_pos0.x, p_force->_pos0.y, p_force->_pos0.z);
					force->lowdim.dir = make_float3 (p_force->_direction.x, p_force->_direction.y, p_force->_direction.z);
					force->lowdim.visX = p_force->_visX;
					force->lowdim.visY = p_force->_visY;
					force->lowdim.visZ = p_force->_visZ;
				}
				else if (typeid (*(p->ext_forces[j])) == typeid (repulsion_plane) ) {
					RepulsionPlane<number> * p_force = (RepulsionPlane<number> *) p->ext_forces[j];
					force->type = CUDA_REPULSION_PLANE;
					force->repulsionplane.stiff = p_force->_stiff;
					force->repulsionplane.position = p_force->_position;
					force->repulsionplane.dir = make_float3 (p_force->_direction.x, p_force->_direction.y, p_force->_direction.z);
				}
				else if (typeid (*(p->ext_forces[j])) == typeid (repulsion_plane_moving) ) {
					RepulsionPlaneMoving<number> * p_force = (RepulsionPlaneMoving<number> *) p->ext_forces[j];
					force->type = CUDA_REPULSION_PLANE_MOVING;
					force->repulsionplanemoving.stiff = p_force->_stiff;
					force->repulsionplanemoving.dir = make_float3 (p_force->_direction.x, p_force->_direction.y, p_force->_direction.z);
					force->repulsionplanemoving.p_ind = p_force->_p_ptr->index;

				}
				else if (typeid (*(p->ext_forces[j])) == typeid (LJ_wall) ) {
					LJWall<number> * p_force = (LJWall<number> *) p->ext_forces[j];
					force->type = CUDA_LJ_WALL;
					force->ljwall.stiff = p_force->_stiff;
					force->ljwall.position = p_force->_position;
					force->ljwall.n = p_force->_n;
					force->ljwall.cutoff = p_force->_cutoff;
					force->ljwall.sigma = p_force->_sigma;
					force->ljwall.dir = make_float3 (p_force->_direction.x, p_force->_direction.y, p_force->_direction.z);

				}
				else {
					throw oxDNAException ("Only ConstantRate, MutualTrap, MovingTrap, LowdimMovingTrap, RepulsionPlane and RepulsionPlaneMoving forces are supported on CUDA at the moment.\n");
				}
			}
		}

		CUDA_SAFE_CALL( cudaMemcpy(_d_ext_forces, _h_ext_forces, this->_N * MAX_EXT_FORCES * sizeof (CUDA_trap<number>), cudaMemcpyHostToDevice) );
	}

	// used in the hilbert curve sorting
	if(this->_sort_every > 0)	{
		CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number4>(&_d_buff_vels, this->_vec_size) );
		CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number4>(&_d_buff_Ls, this->_vec_size) );
	}

	// these values are changed only if the curve sorting is enabled
	_h_gpu_index = new int[this->_N];
	_h_cpu_index = new int[this->_N];
	for(int i = 0; i < this->_N; i++) {
		_h_gpu_index[i] = i;
		_h_cpu_index[i] = i;
	}

	for(int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = this->_particles[i];
		if(p->is_rigid_body()) _any_rigid_body = true;
	}

	// copy all the particle related stuff and the constants to device memory
	_host_particles_to_gpu();
	_init_CUDA_MD_symbols();

	_cuda_thermostat->set_seed(lrand48());
	_cuda_thermostat->init(this->_N);

	OX_DEBUG("Allocated CUDA memory: %.2lf MBs", GpuUtils::get_allocated_mem_mb());

	// initialise lists and compute the forces for the first step
	this->_cuda_lists->update(this->_d_poss, this->_d_list_poss, this->_d_bonds);
	_curr_step = this->_read_conf_step;
	if (_restart_step_counter) _curr_step = 0;
	_set_external_forces();
	this->_cuda_interaction->compute_forces(this->_cuda_lists, this->_d_poss, this->_d_orientations, _d_forces, _d_torques, this->_d_bonds);
}

template<typename number, typename number4>
void MD_CUDABackend<number, number4>::print_conf(llint curr_step, bool reduced, bool only_last) {
	_gpu_to_host_particles();
	MDBackend<number>::print_conf(curr_step, reduced, only_last);
}

template<typename number, typename number4>
void MD_CUDABackend<number, number4>::_print_ready_observables(llint curr_step) {
	_gpu_to_host_particles();
	this->_lists->global_update(true);
	MDBackend<number>::_print_ready_observables(curr_step);
	_host_particles_to_gpu();
}

template<typename number, typename number4>
void MD_CUDABackend<number, number4>::fix_diffusion() {
	_gpu_to_host_particles();
	MDBackend<number>::fix_diffusion();
	_host_particles_to_gpu();
}

// template instantiations
template class MD_CUDABackend<float, float4>;
template class MD_CUDABackend<double, LR_double4>;
