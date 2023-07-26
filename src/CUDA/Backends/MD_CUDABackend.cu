/*
 * MD_CUDABackend.cpp
 *
 *  Created on: 03/set/2010
 *      Author: lorenzo
 */

#include "MD_CUDABackend.h"

#include "../CUDAForces.h"
#include "CUDA_MD.cuh"
#include "../CUDA_base_interactions.h"
#include "../../Interactions/DNAInteraction.h"
#include "../../Observables/ObservableOutput.h"
#include "../Thermostats/CUDAThermostatFactory.h"

#include <typeinfo>

// these pragma instructions remove a few nvcc warnings
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvla"

MD_CUDABackend::MD_CUDABackend() :
				MDBackend(),
				CUDABaseBackend(),
				_max_ext_forces(0),
				_error_conf_file("error_conf.dat") {
	_use_edge = false;
	_any_rigid_body = false;

	_d_vels = _d_Ls = _d_forces = _d_torques = _d_buff_vels = nullptr;
	_h_vels = _h_Ls = _h_forces = _h_torques = _d_buff_Ls = nullptr;
	_h_gpu_index = _h_cpu_index = nullptr;

	_d_particles_to_mols = _d_mol_sizes = nullptr;
	_d_molecular_coms = nullptr;
	_d_buff_particles_to_mols = nullptr;

	_d_ext_forces = nullptr;

	_avoid_cpu_calculations = false;

	_timer_sorting = nullptr;

	_barostat_attempts = _barostat_accepted = 0;

	_print_energy = false;

	_obs_output_error_conf = nullptr;

	// on CUDA the timers need to be told to explicitly synchronise on the GPU
	TimingManager::instance()->enable_sync();
}

MD_CUDABackend::~MD_CUDABackend() {
	if(_d_particles_to_mols != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_particles_to_mols));
	}

	if(_d_vels != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_vels));
		CUDA_SAFE_CALL(cudaFree(_d_Ls));
		CUDA_SAFE_CALL(cudaFree(_d_forces));
		CUDA_SAFE_CALL(cudaFree(_d_torques));
	}

	if(_d_molecular_coms != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_molecular_coms));
	}

	if(_sort_every > 0 && _d_buff_vels != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_buff_vels));
		CUDA_SAFE_CALL(cudaFree(_d_buff_Ls));
		CUDA_SAFE_CALL(cudaFree(_d_buff_particles_to_mols));
	}

	if(_h_gpu_index != nullptr) {
		delete[] _h_gpu_index;
		delete[] _h_cpu_index;
	}

	if(_external_forces) {
		if(_d_ext_forces != nullptr) {
			CUDA_SAFE_CALL(cudaFree(_d_ext_forces));
		}
	}

	if(_h_vels != nullptr) {
		delete[] _h_vels;
		delete[] _h_Ls;
		delete[] _h_forces;
		delete[] _h_torques;
	}

	if(_obs_output_error_conf != nullptr) {
		delete _obs_output_error_conf;
	}
}

void MD_CUDABackend::_host_to_gpu() {
	CUDABaseBackend::_host_to_gpu();
	CUDA_SAFE_CALL(cudaMemcpy(_d_particles_to_mols, _h_particles_to_mols.data(), sizeof(int) * N(), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(_d_vels, _h_vels, _vec_size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(_d_Ls, _h_Ls, _vec_size, cudaMemcpyHostToDevice));
}

void MD_CUDABackend::_gpu_to_host() {
	CUDABaseBackend::_gpu_to_host();
	CUDA_SAFE_CALL(cudaMemcpy(_h_particles_to_mols.data(), _d_particles_to_mols, sizeof(int) * N(), cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpy(_h_vels, _d_vels, _vec_size, cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpy(_h_Ls, _d_Ls, _vec_size, cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpy(_h_forces, _d_forces, _vec_size, cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpy(_h_torques, _d_torques, _vec_size, cudaMemcpyDeviceToHost));
}

void MD_CUDABackend::_apply_external_forces_changes() {
	if(_external_forces) {
		if(_sort_every > 0) {
			throw oxDNAException("External forces and CUDA_sort_every > 0 are not compatible");
		}
		bool first_time = (_d_ext_forces == nullptr);

		static std::vector<CUDA_trap> h_ext_forces(N() * MAX_EXT_FORCES);

		if(first_time) {
			// TODO: possible memory leak on oxpy
			CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<CUDA_trap >(&_d_ext_forces, N() * MAX_EXT_FORCES * sizeof(CUDA_trap)));

			for(int i = 0; i < N() * MAX_EXT_FORCES; i++) {
				h_ext_forces[i].type = -1;
			}
		}

		for(int i = 0; i < N(); i++) {
			BaseParticle *p = _particles[i];

			for(uint j = 0; j < p->ext_forces.size(); j++) {
				_max_ext_forces = max(_max_ext_forces, (int) p->ext_forces.size());
				auto &force_type = typeid(*(p->ext_forces[j]));

				CUDA_trap *cuda_force = &(h_ext_forces[j * N() + i]);
				if(force_type == typeid(ConstantRateForce)) {
					ConstantRateForce *p_force = (ConstantRateForce *) p->ext_forces[j];
					init_ConstantRateForce_from_CPU(&cuda_force->constant, p_force);
				}
				else if(force_type == typeid(MutualTrap)) {
					MutualTrap *p_force = (MutualTrap *) p->ext_forces[j];
					init_MutualTrap_from_CPU(&cuda_force->mutual, p_force);
				}
				else if(force_type == typeid(MovingTrap)) {
					MovingTrap *p_force = (MovingTrap *) p->ext_forces[j];
					init_MovingTrap_from_CPU(&cuda_force->moving, p_force);
				}
				else if(force_type == typeid(LowdimMovingTrap)) {
					LowdimMovingTrap *p_force = (LowdimMovingTrap *) p->ext_forces[j];
					init_LowdimMovingTrap_from_CPU(&cuda_force->lowdim, p_force);
				}
				else if(force_type == typeid(RepulsionPlane)) {
					RepulsionPlane *p_force = (RepulsionPlane *) p->ext_forces[j];
					init_RepulsionPlane_from_CPU(&cuda_force->repulsionplane, p_force);
				}
				else if(force_type == typeid(RepulsionPlaneMoving)) {
					RepulsionPlaneMoving *p_force = (RepulsionPlaneMoving *) p->ext_forces[j];
					init_RepulsionPlaneMoving_from_CPU(&cuda_force->repulsionplanemoving, p_force);
				}
				else if(force_type == typeid(RepulsiveSphere)) {
					RepulsiveSphere *p_force = (RepulsiveSphere *) p->ext_forces[j];
					init_RepulsiveSphere_from_CPU(&cuda_force->repulsivesphere, p_force);
				}
				else if(force_type == typeid(RepulsiveSphereSmooth)) {
					RepulsiveSphereSmooth *p_force = (RepulsiveSphereSmooth *) p->ext_forces[j];
					init_RepulsiveSphereSmooth_from_CPU(&cuda_force->repulsivespheresmooth, p_force);
				}
				else if(force_type == typeid(LJWall)) {
					LJWall *p_force = (LJWall *) p->ext_forces[j];
					init_LJWall_from_CPU(&cuda_force->ljwall, p_force);
				}
				else if(force_type == typeid(ConstantRateTorque)) {
					ConstantRateTorque *p_force = (ConstantRateTorque *) p->ext_forces[j];
					init_ConstantRateTorque_from_CPU(&cuda_force->constantratetorque, p_force);
				}
				else if(force_type == typeid(GenericCentralForce)) {
					GenericCentralForce *p_force = (GenericCentralForce *) p->ext_forces[j];
					init_GenericCentralForce_from_CPU(&cuda_force->genericconstantforce, p_force);
				}
				else if(force_type == typeid(LJCone)) {
					LJCone *p_force = (LJCone *) p->ext_forces[j];
					init_LJCone_from_CPU(&cuda_force->ljcone, p_force);
				}
				else if(force_type == typeid(RepulsiveEllipsoid)) {
					RepulsiveEllipsoid *p_force = (RepulsiveEllipsoid *) p->ext_forces[j];
					init_RepulsiveEllipsoid_from_CPU(&cuda_force->repulsiveellipsoid, p_force);
				}
				else if(force_type == typeid(COMForce)) {
					COMForce *p_force = (COMForce *) p->ext_forces[j];
					init_COMForce_from_CPU(&cuda_force->comforce, p_force, first_time);
				}
				else if(force_type == typeid(LTCOMTrap)) {
					LTCOMTrap *p_force = (LTCOMTrap *) p->ext_forces[j];
					init_LTCOMTrap_from_CPU(&cuda_force->ltcomtrap, p_force, first_time);
				}
				else {
					throw oxDNAException("Only ConstantRate, MutualTrap, MovingTrap, LowdimMovingTrap, RepulsionPlane, "
							"RepulsionPlaneMoving, RepulsiveSphere, LJWall, ConstantRateTorque, GenericCentralForce, "
							"RepulsiveEllipsoid, COMForce and LTCOMTrap"
							"forces are supported on CUDA at the moment.\n");
				}
			}
		}

		CUDA_SAFE_CALL(cudaMemcpy(_d_ext_forces, h_ext_forces.data(), N() * MAX_EXT_FORCES * sizeof(CUDA_trap), cudaMemcpyHostToDevice));
	}
}

void MD_CUDABackend::apply_changes_to_simulation_data() {
	for(int i = 0; i < N(); i++) {
		int gpu_index = _h_gpu_index[i];
		BaseParticle *p = _particles[gpu_index];

		_h_poss[i].x = p->pos.x;
		_h_poss[i].y = p->pos.y;
		_h_poss[i].z = p->pos.z;

		_h_particles_to_mols[i] = p->strand_id;

		// convert index and type into a float
		int msk = -1 << 22; // binary mask: all 1's and 22 0's;
		// btype has a sign, and thus has to go first
		_h_poss[i].w = GpuUtils::int_as_float((p->btype << 22) | ((~msk) & p->index));
		// we immediately check that the index and base type that we read are sensible
		int mybtype = (GpuUtils::float_as_int(_h_poss[i].w)) >> 22;
		int myindex = (GpuUtils::float_as_int(_h_poss[i].w)) & (~msk);
		if(p->btype != mybtype) {
			throw oxDNAException("Could not treat the type (A, C, G, T or something specific) of particle %d; On CUDA, integer base types cannot be larger than 511 or smaller than -511");
		}
		if(p->index != myindex) {
			throw oxDNAException("Could not treat the index of particle %d; remember that on CUDA the maximum c_number of particles is 2^21", p->index);
		}
		
		if(p->n3 == P_VIRTUAL) {
			_h_bonds[i].n3 = P_INVALID;
		}
		else {
			_h_bonds[i].n3 = _h_cpu_index[p->n3->index];
		}
		if(p->n5 == P_VIRTUAL) {
			_h_bonds[i].n5 = P_INVALID;
		}
		else {
			_h_bonds[i].n5 = _h_cpu_index[p->n5->index];
		}

		_h_vels[i].x = p->vel.x;
		_h_vels[i].y = p->vel.y;
		_h_vels[i].z = p->vel.z;

		_h_Ls[i].x = p->L.x;
		_h_Ls[i].y = p->L.y;
		_h_Ls[i].z = p->L.z;

		number trace = p->orientation.v1.x + p->orientation.v2.y + p->orientation.v3.z;
		if(trace > 0) {
			number s = 0.5 / sqrt(trace + 1.0);
			_h_orientations[i].w = 0.25 / s;
			_h_orientations[i].x = (p->orientation.v3.y - p->orientation.v2.z) * s;
			_h_orientations[i].y = (p->orientation.v1.z - p->orientation.v3.x) * s;
			_h_orientations[i].z = (p->orientation.v2.x - p->orientation.v1.y) * s;
		}
		else { // Finding largest diagonal element
			if((p->orientation.v1.x > p->orientation.v2.y) && (p->orientation.v1.x > p->orientation.v3.z)) {
				number s = 0.5 / sqrt(1.0 + p->orientation.v1.x - p->orientation.v2.y - p->orientation.v3.z);
				_h_orientations[i].w = (p->orientation.v3.y - p->orientation.v2.z) * s;
				_h_orientations[i].x = 0.25 / s;
				_h_orientations[i].y = (p->orientation.v1.y + p->orientation.v2.x) * s;
				_h_orientations[i].z = (p->orientation.v1.z + p->orientation.v3.x) * s;
			}
			else if(p->orientation.v2.y > p->orientation.v3.z) {
				number s = 0.5 / sqrt(1.0 + p->orientation.v2.y - p->orientation.v1.x - p->orientation.v3.z);
				_h_orientations[i].w = (p->orientation.v1.z - p->orientation.v3.x) * s;
				_h_orientations[i].x = (p->orientation.v1.y + p->orientation.v2.x) * s;
				_h_orientations[i].y = 0.25 / s;
				_h_orientations[i].z = (p->orientation.v2.z + p->orientation.v3.y) * s;
			}
			else {
				number s = 0.5 / sqrt(1.0 + p->orientation.v3.z - p->orientation.v1.x - p->orientation.v2.y);
				_h_orientations[i].w = (p->orientation.v2.x - p->orientation.v1.y) * s;
				_h_orientations[i].x = (p->orientation.v1.z + p->orientation.v3.x) * s;
				_h_orientations[i].y = (p->orientation.v2.z + p->orientation.v3.y) * s;
				_h_orientations[i].z = 0.25 / s;
			}
		}
	}
	_host_to_gpu();
	_apply_external_forces_changes();

	_cuda_interaction->sync_GPU();
}

void MD_CUDABackend::apply_simulation_data_changes() {
	_gpu_to_host();

	for(int i = 0; i < N(); i++) {
		// since we may have been sorted all the particles in a different order
		// we first take the particle index from the 4th component of its
		// position, and then use that index to access the right BaseParticle pointer
		int msk = (-1 << 22);
		int newindex = ((GpuUtils::float_as_int(_h_poss[i].w)) & (~msk));
		_h_gpu_index[i] = newindex;
		_h_cpu_index[newindex] = i;
		BaseParticle *p = _particles[newindex];
		assert(p->index == newindex);

		p->pos.x = _h_poss[i].x;
		p->pos.y = _h_poss[i].y;
		p->pos.z = _h_poss[i].z;

		p->strand_id = _h_particles_to_mols[i];

		// get index and type from the fourth component of the position
		p->btype = (GpuUtils::float_as_int(_h_poss[i].w)) >> 22;

		if(_h_bonds[i].n3 == P_INVALID) {
			p->n3 = P_VIRTUAL;
		}
		else {
			int n3index = ((GpuUtils::float_as_int(_h_poss[_h_bonds[i].n3].w)) & (~msk));
			p->n3 = _particles[n3index];
		}
		if(_h_bonds[i].n5 == P_INVALID) {
			p->n5 = P_VIRTUAL;
		}
		else {
			int n5index = ((GpuUtils::float_as_int(_h_poss[_h_bonds[i].n5].w)) & (~msk));
			p->n5 = _particles[n5index];
		}

		p->vel.x = _h_vels[i].x;
		p->vel.y = _h_vels[i].y;
		p->vel.z = _h_vels[i].z;

		p->L.x = _h_Ls[i].x;
		p->L.y = _h_Ls[i].y;
		p->L.z = _h_Ls[i].z;

		// we cast c_numbers (which are usually floats) to numbers (which are usually double) to increase the numerical accuracy
		number sqx = (number) _h_orientations[i].x * (number) _h_orientations[i].x;
		number sqy = (number) _h_orientations[i].y * (number) _h_orientations[i].y;
		number sqz = (number) _h_orientations[i].z * (number) _h_orientations[i].z;
		number sqw = (number) _h_orientations[i].w * (number) _h_orientations[i].w;
		number xy = (number) _h_orientations[i].x * (number) _h_orientations[i].y;
		number xz = (number) _h_orientations[i].x * (number) _h_orientations[i].z;
		number xw = (number) _h_orientations[i].x * (number) _h_orientations[i].w;
		number yz = (number) _h_orientations[i].y * (number) _h_orientations[i].z;
		number yw = (number) _h_orientations[i].y * (number) _h_orientations[i].w;
		number zw = (number) _h_orientations[i].z * (number) _h_orientations[i].w;
		number invs = 1.0 / (sqx + sqy + sqz + sqw);

		p->orientation.v1.x = (sqx - sqy - sqz + sqw) * invs;
		p->orientation.v1.y = 2 * (xy - zw) * invs;
		p->orientation.v1.z = 2 * (xz + yw) * invs;
		p->orientation.v2.x = 2 * (xy + zw) * invs;
		p->orientation.v2.y = (-sqx + sqy - sqz + sqw) * invs;
		p->orientation.v2.z = 2 * (yz - xw) * invs;
		p->orientation.v3.x = 2 * (xz - yw) * invs;
		p->orientation.v3.y = 2 * (yz + xw) * invs;
		p->orientation.v3.z = (-sqx - sqy + sqz + sqw) * invs;

		p->orientationT = p->orientation.get_transpose();
		p->set_positions();
	}

	_cuda_interaction->sync_host();

	if(!_avoid_cpu_calculations) {
		_lists->global_update(true);
	}
}

void MD_CUDABackend::_init_CUDA_MD_symbols() {
	float f_copy = _sqr_verlet_skin;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_sqr_verlet_skin, &f_copy, sizeof(float)));
	f_copy = _dt;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_dt, &f_copy, sizeof(float)));
	int myN = N();
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &myN, sizeof(int)));
}

void MD_CUDABackend::_first_step() {
	first_step
		<<<_particles_kernel_cfg.blocks, _particles_kernel_cfg.threads_per_block>>>
		(_d_poss, _d_orientations, _d_list_poss, _d_vels, _d_Ls, _d_forces, _d_torques, _d_are_lists_old);
	CUT_CHECK_ERROR("_first_step error");
}

void MD_CUDABackend::_rescale_molecular_positions(c_number4 new_Ls, c_number4 old_Ls, bool recompute_coms) {
	c_number4 shift_factor = {
		new_Ls.x / old_Ls.x - (c_number) 1.f,
		new_Ls.y / old_Ls.y - (c_number) 1.f,
		new_Ls.z / old_Ls.z - (c_number) 1.f,
		0.
	};

	if(recompute_coms) {
		CUDA_SAFE_CALL(cudaMemset(_d_molecular_coms, 0, sizeof(c_number4) * _molecules.size()));
		compute_molecular_coms
			<<<_particles_kernel_cfg.blocks, _particles_kernel_cfg.threads_per_block>>>
			(_d_molecular_coms, _d_particles_to_mols, _d_mol_sizes, _d_poss);
		CUT_CHECK_ERROR("compute_molecular_coms");
	}
	else {
		shift_factor = -shift_factor;
	}

	rescale_molecular_positions
		<<<_particles_kernel_cfg.blocks, _particles_kernel_cfg.threads_per_block>>>
		(_d_molecular_coms, _d_particles_to_mols, _d_poss, shift_factor);
	CUT_CHECK_ERROR("rescale_molecular_positions");
}

void MD_CUDABackend::_rescale_positions(c_number4 new_Ls, c_number4 old_Ls) {
	c_number4 ratio = {
		new_Ls.x / old_Ls.x,
		new_Ls.y / old_Ls.y,
		new_Ls.z / old_Ls.z,
		0.
	};

	rescale_positions
		<<<_particles_kernel_cfg.blocks, _particles_kernel_cfg.threads_per_block>>>
		(_d_poss, ratio);
	CUT_CHECK_ERROR("_rescale_positions error");
}

void MD_CUDABackend::_apply_barostat() {
	_barostat_attempts++;

	_set_external_forces();
	_cuda_interaction->compute_forces(_cuda_lists, _d_poss, _d_orientations, _d_forces, _d_torques, _d_bonds, _d_cuda_box);
	double old_energy = GpuUtils::sum_c_number4_to_double_on_GPU(_d_forces, N()) / 2.;
	c_number old_V = _h_cuda_box.V();
	c_number4 old_Ls = _h_cuda_box.box_sides();

	c_number4 new_Ls = old_Ls;
	if(_barostat_isotropic) {
		c_number dL = _delta_L * (drand48() - (c_number) 0.5);
		new_Ls.x += dL;
		new_Ls.y += dL;
		new_Ls.z += dL;
	}
	else {
		new_Ls.x += _delta_L * (drand48() - (c_number) 0.5);
		new_Ls.y += _delta_L * (drand48() - (c_number) 0.5);
		new_Ls.z += _delta_L * (drand48() - (c_number) 0.5);
	}
	_h_cuda_box.change_sides(new_Ls.x, new_Ls.y, new_Ls.z);
	CUDA_SAFE_CALL(cudaMemcpy(_d_cuda_box, &_h_cuda_box, sizeof(CUDABox), cudaMemcpyHostToDevice));

	if(_barostat_molecular) {
		_rescale_molecular_positions(new_Ls, old_Ls, true);
	}
	else {
		_rescale_positions(new_Ls, old_Ls);
	}
	_cuda_lists->update(_d_poss, _d_list_poss, _d_bonds);

	_set_external_forces();
	_cuda_interaction->compute_forces(_cuda_lists, _d_poss, _d_orientations, _d_forces, _d_torques, _d_bonds, _d_cuda_box);
	double new_energy = GpuUtils::sum_c_number4_to_double_on_GPU(_d_forces, N()) / 2.;
	c_number new_V = _h_cuda_box.V();

	// acceptance
	c_number dE = new_energy - old_energy;
	c_number dV = new_V - old_V;
	int N_objs = (_barostat_molecular) ? _molecules.size() : N();
	c_number acc = exp(-(dE + _P * dV - N_objs * _T * log(new_V / old_V)) / _T);
	// accepted
	if(acc > drand48()) {
		_barostat_accepted++;
		CONFIG_INFO->notify(CONFIG_INFO->box->UPDATE_EVENT);
	}
	// rejected
	else {
		_h_cuda_box.change_sides(old_Ls.x, old_Ls.y, old_Ls.z);
		CUDA_SAFE_CALL(cudaMemcpy(_d_cuda_box, &_h_cuda_box, sizeof(CUDABox), cudaMemcpyHostToDevice));
		if(_barostat_molecular) {
			_rescale_molecular_positions(new_Ls, old_Ls, false);
		}
		else {
			_rescale_positions(old_Ls, new_Ls);
		}
		_cuda_lists->update(_d_poss, _d_list_poss, _d_bonds);
	}
	_barostat_acceptance = _barostat_accepted / (c_number) _barostat_attempts;

	if(_cuda_barostat_always_refresh) {
		// if the user wishes so, we refresh all the velocities after each barostat attempt
		_cuda_barostat_thermostat->apply_cuda(_d_poss, _d_orientations, _d_vels, _d_Ls, current_step());
	}
}

void MD_CUDABackend::_forces_second_step() {
	_set_external_forces();
	_cuda_interaction->compute_forces(_cuda_lists, _d_poss, _d_orientations, _d_forces, _d_torques, _d_bonds, _d_cuda_box);

	second_step
		<<<_particles_kernel_cfg.blocks, _particles_kernel_cfg.threads_per_block>>>
		(_d_vels, _d_Ls, _d_forces, _d_torques);
	CUT_CHECK_ERROR("second_step");
}

void MD_CUDABackend::_set_external_forces() {
	set_external_forces
		<<<_particles_kernel_cfg.blocks, _particles_kernel_cfg.threads_per_block>>>
		(_d_poss, _d_orientations, _d_ext_forces, _d_forces, _d_torques, current_step(), _max_ext_forces, _d_cuda_box);
	CUT_CHECK_ERROR("set_external_forces");
}

void MD_CUDABackend::_sort_particles() {
	CUDABaseBackend::_sort_index();
	permute_particles
		<<<_particles_kernel_cfg.blocks, _particles_kernel_cfg.threads_per_block>>>
		(_d_sorted_hindex, _d_inv_sorted_hindex, _d_poss, _d_vels, _d_Ls, _d_orientations, _d_bonds, _d_particles_to_mols,
		 _d_buff_poss, _d_buff_vels, _d_buff_Ls, _d_buff_orientations, _d_buff_bonds, _d_buff_particles_to_mols);
	CUT_CHECK_ERROR("_permute_particles error");
	CUDA_SAFE_CALL(cudaMemcpy(_d_orientations, _d_buff_orientations, _orient_size, cudaMemcpyDeviceToDevice));

	// copy back the sorted vectors
	CUDA_SAFE_CALL(cudaMemcpy(_d_poss, _d_buff_poss, _vec_size, cudaMemcpyDeviceToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(_d_bonds, _d_buff_bonds, _bonds_size, cudaMemcpyDeviceToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(_d_vels, _d_buff_vels, _vec_size, cudaMemcpyDeviceToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(_d_Ls, _d_buff_Ls, _vec_size, cudaMemcpyDeviceToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(_d_particles_to_mols, _d_buff_particles_to_mols, sizeof(int) * N(), cudaMemcpyDeviceToDevice));
}

void MD_CUDABackend::_thermalize() {
	_cuda_thermostat->apply_cuda(_d_poss, _d_orientations, _d_vels, _d_Ls, current_step());
}

void MD_CUDABackend::_update_stress_tensor() {
	if(_update_st_every > 0 && (CONFIG_INFO->curr_step % _update_st_every == 0)) {
		_interaction->set_stress_tensor(_cuda_interaction->CPU_stress_tensor(_d_vels));
	}
}

void MD_CUDABackend::sim_step() {
	_mytimer->resume();

	_timer_first_step->resume();
	_first_step();
	_timer_first_step->pause();

	_timer_sorting->resume();
	if(_d_are_lists_old[0] && _sort_every > 0 && (_N_updates % _sort_every == 0)) {
		_sort_particles();
	}
	_timer_sorting->pause();

	_timer_lists->resume();
	if(_d_are_lists_old[0]) {
		try {
			_cuda_lists->update(_d_poss, _d_list_poss, _d_bonds);
		}
		catch (oxDNAException &e) {
			apply_simulation_data_changes();
			_obs_output_error_conf->print_output(current_step());
			throw oxDNAException("%s ----> The last configuration has been printed to %s", e.what(), _error_conf_file.c_str());
		}
		_d_are_lists_old[0] = false;
		_N_updates++;
	}
	_timer_lists->pause();

	if(_is_barostat_active()) {
		_timer_barostat->resume();
		_apply_barostat();
		_timer_barostat->pause();
	}

	_timer_forces->resume();
	_forces_second_step();
	if(_print_energy) {
		c_number energy = GpuUtils::sum_c_number4_to_double_on_GPU(_d_forces, N());
		_backend_info = Utils::sformat("\tCUDA_energy: %lf", energy / (2. * N()));
	}

	_update_stress_tensor();

	_timer_forces->pause();

	_timer_thermostat->resume();
	_thermalize();
	_timer_thermostat->pause();

	_mytimer->pause();
}

void MD_CUDABackend::get_settings(input_file &inp) {
	MDBackend::get_settings(inp);
	CUDABaseBackend::get_settings(inp);

	if(getInputBool(&inp, "use_edge", &_use_edge, 0) == KEY_FOUND) {
		if(_use_edge && sizeof(c_number) == sizeof(double)) {
			throw oxDNAException("use_edge and double precision are not compatible");
		}
		if(_use_edge && _use_barostat) {
			throw oxDNAException("use_edge and use_barostat are not compatible");
		}
	}

	getInputBool(&inp, "CUDA_avoid_cpu_calculations", &_avoid_cpu_calculations, 0);
	getInputBool(&inp, "CUDA_barostat_always_refresh", &_cuda_barostat_always_refresh, 0);
	getInputBool(&inp, "CUDA_print_energy", &_print_energy, 0);
	getInputInt(&inp, "CUDA_update_stress_tensor_every", &_update_st_every, 0);

	_cuda_thermostat = CUDAThermostatFactory::make_thermostat(inp, _box.get());
	_cuda_thermostat->get_settings(inp);

	if(_use_barostat) {
		std::string input_string = Utils::sformat("newtonian_steps = 1\npt = 1.0\ndt = %lf\nT = %lf", this->_dt, this->_T);
		input_file *inp_file = Utils::get_input_file_from_string(input_string);
		_cuda_barostat_thermostat = std::make_shared<CUDABrownianThermostat>();
		_cuda_barostat_thermostat->get_settings(*inp_file);
	}

	std::string init_string = Utils::sformat("{\n\tname = %s\n\tprint_every = 0\n\tonly_last = 1\n}\n", _error_conf_file.c_str());
	_obs_output_error_conf = new ObservableOutput(init_string);
	_obs_output_error_conf->add_observable("type = configuration");

	// if we want to limit the calculations done on CPU we clear the default ObservableOutputs and tell them to just print the timesteps (and, for constant-pressure, simulations, also the density)
	if(_avoid_cpu_calculations) {
		_obs_output_file->clear();
		_obs_output_file->add_observable("type = step\nunits = MD");
		if(_use_barostat) {
			_obs_output_file->add_observable("type = density");
		}

		bool no_stdout_energy = false;
		getInputBool(&inp, "no_stdout_energy", &no_stdout_energy, 0);
		if(!no_stdout_energy) {
			_obs_output_stdout->clear();
			_obs_output_stdout->add_observable("type = step");
			_obs_output_stdout->add_observable("type = step\nunits = MD");
			if(_use_barostat) {
				_obs_output_stdout->add_observable("type = density");
			}
		}
	}
}

void MD_CUDABackend::init() {
	MDBackend::init();
	CUDABaseBackend::init_cuda();

	_timer_sorting = TimingManager::instance()->new_timer(std::string("Hilbert sorting"), std::string("SimBackend"));

	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<c_number4>(&_d_vels, _vec_size));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<c_number4>(&_d_Ls, _vec_size));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<c_number4>(&_d_forces, _vec_size));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<c_number4>(&_d_torques, _vec_size));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&_d_particles_to_mols, sizeof(int) * N()));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&_d_mol_sizes, sizeof(int) * _molecules.size()));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<c_number4>(&_d_molecular_coms, sizeof(c_number4) * _molecules.size()));

	CUDA_SAFE_CALL(cudaMemset(_d_forces, 0, _vec_size));
	CUDA_SAFE_CALL(cudaMemset(_d_torques, 0, _vec_size));

	_h_particles_to_mols.resize(N());
	_h_vels = new c_number4[N()];
	_h_Ls = new c_number4[N()];
	_h_forces = new c_number4[N()];
	_h_torques = new c_number4[N()];

	_obs_output_error_conf->init();

	// initialise the GPU array containing the size of the molecules
	std::vector<int> mol_sizes;
	for(auto mol : _molecules) {
		mol_sizes.push_back(mol->N());
	}
	CUDA_SAFE_CALL(cudaMemcpy(_d_mol_sizes, mol_sizes.data(), sizeof(int) * _molecules.size(), cudaMemcpyHostToDevice));

	_apply_external_forces_changes();

	// used in the hilbert curve sorting
	if(_sort_every > 0) {
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<c_number4>(&_d_buff_vels, _vec_size));
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<c_number4>(&_d_buff_Ls, _vec_size));
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&_d_buff_particles_to_mols, sizeof(int) * N()));
	}

	// these values are changed only if the curve sorting is enabled
	_h_gpu_index = new int[N()];
	_h_cpu_index = new int[N()];
	for(int i = 0; i < N(); i++) {
		_h_gpu_index[i] = i;
		_h_cpu_index[i] = i;
	}

	for(int i = 0; i < N(); i++) {
		BaseParticle *p = _particles[i];
		if(p->is_rigid_body()) {
			_any_rigid_body = true;
		}
	}

	// copy all the particle related stuff and the constants to device memory
	apply_changes_to_simulation_data();
	_init_CUDA_MD_symbols();

	_cuda_thermostat->set_seed(lrand48());
	_cuda_thermostat->init();

	if(_use_barostat) {
		_cuda_barostat_thermostat->set_seed(lrand48());
		_cuda_barostat_thermostat->init();
	}

	OX_LOG(Logger::LOG_INFO, "Allocated CUDA memory: %.2lf MBs", GpuUtils::get_allocated_mem_mb());

	// initialise lists and compute the forces for the first step
	_cuda_lists->update(_d_poss, _d_list_poss, _d_bonds);
	_set_external_forces();
	_cuda_interaction->compute_forces(_cuda_lists, _d_poss, _d_orientations, _d_forces, _d_torques, _d_bonds, _d_cuda_box);

	if(_update_st_every > 0) {
		_interaction->set_stress_tensor(_cuda_interaction->CPU_stress_tensor(_d_vels));
	}
}

#pragma GCC diagnostic pop
