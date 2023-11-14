/*
 * CUDAForces.h
 *
 *  Created on: 26/giu/2013
 *      Author: lorenzo
 */

#ifndef CUDAFORCES_H_
#define CUDAFORCES_H_

#include "../Forces/COMForce.h"
#include "../Forces/ConstantRateForce.h"
#include "../Forces/ConstantRateTorque.h"
#include "../Forces/ConstantTrap.h"
#include "../Forces/LowdimMovingTrap.h"
#include "../Forces/MovingTrap.h"
#include "../Forces/MutualTrap.h"
#include "../Forces/RepulsionPlane.h"
#include "../Forces/RepulsionPlaneMoving.h"
#include "../Forces/RepulsiveSphere.h"
#include "../Forces/RepulsiveSphereSmooth.h"
#include "../Forces/LJWall.h"
#include "../Forces/GenericCentralForce.h"
#include "../Forces/LJCone.h"
#include "../Forces/RepulsiveEllipsoid.h"
#include "../Forces/Metadynamics/LTCOMTrap.h"

#include "CUDAUtils.h"

#define CUDA_TRAP_NO_FORCE -1
#define CUDA_TRAP_CONSTANT 0
#define CUDA_TRAP_MUTUAL 1
#define CUDA_TRAP_MOVING 2
#define CUDA_REPULSION_PLANE 3
#define CUDA_REPULSION_PLANE_MOVING 4
#define CUDA_TRAP_MOVING_LOWDIM 5
#define CUDA_LJ_WALL 6
#define CUDA_REPULSIVE_SPHERE 7
#define CUDA_CONSTANT_RATE_TORQUE 8
#define CUDA_GENERIC_CENTRAL_FORCE 9
#define CUDA_LJ_CONE 10
#define CUDA_REPULSIVE_SPHERE_SMOOTH 11
#define CUDA_REPULSIVE_ELLIPSOID 12
#define CUDA_COM_FORCE 13
#define CUDA_LR_COM_TRAP 14

/**
 * @brief CUDA version of a ConstantRateForce.
 */
struct constant_rate_force {
	int type;
	c_number x, y, z;
	bool dir_as_centre;
	c_number rate;
	c_number F0;
};

void init_ConstantRateForce_from_CPU(constant_rate_force *cuda_force, ConstantRateForce *cpu_force) {
	cuda_force->type = CUDA_TRAP_CONSTANT;
	cuda_force->F0 = cpu_force->_F0;
	cuda_force->dir_as_centre = cpu_force->dir_as_centre;
	cuda_force->rate = cpu_force->_rate;
	cuda_force->x = cpu_force->_direction.x;
	cuda_force->y = cpu_force->_direction.y;
	cuda_force->z = cpu_force->_direction.z;
}

/**
 * @brief CUDA version of a MutualTrap.
 */
struct mutual_trap {
	int type;
	c_number stiff;
	c_number r0;
	int p_ind;
	bool PBC;
	c_number rate;
	c_number stiff_rate;
};

void init_MutualTrap_from_CPU(mutual_trap *cuda_force, MutualTrap *cpu_force) {
	cuda_force->type = CUDA_TRAP_MUTUAL;
	cuda_force->rate = cpu_force->_rate;
	cuda_force->stiff = cpu_force->_stiff;
	cuda_force->stiff_rate = cpu_force->_stiff_rate;
	cuda_force->r0 = cpu_force->_r0;
	cuda_force->p_ind = cpu_force->_p_ptr->index;
	cuda_force->PBC = cpu_force->PBC;
}

/**
 * @brief CUDA version of a MovingTrap.
 */
struct moving_trap {
	int type;
	c_number stiff;
	c_number rate;
	float3 pos0;
	float3 dir;
};

void init_MovingTrap_from_CPU(moving_trap *cuda_force, MovingTrap *cpu_force) {
	cuda_force->type = CUDA_TRAP_MOVING;
	cuda_force->stiff = cpu_force->_stiff;
	cuda_force->rate = cpu_force->_rate;
	cuda_force->pos0 = make_float3(cpu_force->_pos0.x, cpu_force->_pos0.y, cpu_force->_pos0.z);
	cuda_force->dir = make_float3(cpu_force->_direction.x, cpu_force->_direction.y, cpu_force->_direction.z);
}

/**
 * @brief CUDA version of a LowdimMovingTrap.
 */
struct lowdim_moving_trap {
	int type;
	c_number stiff;
	c_number rate;
	float3 pos0;
	float3 dir;
	bool visX;
	bool visY;
	bool visZ;
};

void init_LowdimMovingTrap_from_CPU(lowdim_moving_trap *cuda_force, LowdimMovingTrap *cpu_force) {
	cuda_force->type = CUDA_TRAP_MOVING_LOWDIM;
	cuda_force->stiff = cpu_force->_stiff;
	cuda_force->rate = cpu_force->_rate;
	cuda_force->pos0 = make_float3(cpu_force->_pos0.x, cpu_force->_pos0.y, cpu_force->_pos0.z);
	cuda_force->dir = make_float3(cpu_force->_direction.x, cpu_force->_direction.y, cpu_force->_direction.z);
	cuda_force->visX = cpu_force->_visX;
	cuda_force->visY = cpu_force->_visY;
	cuda_force->visZ = cpu_force->_visZ;
}

/**
 * @brief CUDA version of a RepulsionPlane.
 */
struct repulsion_plane {
	int type;
	c_number stiff;
	c_number position;
	float3 dir;
};

void init_RepulsionPlane_from_CPU(repulsion_plane *cuda_force, RepulsionPlane *cpu_force) {
	cuda_force->type = CUDA_REPULSION_PLANE;
	cuda_force->stiff = cpu_force->_stiff;
	cuda_force->position = cpu_force->_position;
	cuda_force->dir = make_float3(cpu_force->_direction.x, cpu_force->_direction.y, cpu_force->_direction.z);
}

/**
 * @brief CUDA version of a RepulsionPlaneMoving.
 */
struct repulsion_plane_moving {
	int type;
	c_number stiff;
	int low_idx, high_idx;
	float3 dir;
};

void init_RepulsionPlaneMoving_from_CPU(repulsion_plane_moving *cuda_force, RepulsionPlaneMoving *cpu_force) {
	cuda_force->type = CUDA_REPULSION_PLANE_MOVING;
	cuda_force->stiff = cpu_force->_stiff;
	cuda_force->dir = make_float3(cpu_force->_direction.x, cpu_force->_direction.y, cpu_force->_direction.z);
	cuda_force->low_idx = cpu_force->low_idx;
	cuda_force->high_idx = cpu_force->high_idx;
}

/**
 * @brief CUDA version of a RepulsiveSphere.
 */
struct repulsive_sphere {
	int type;
	c_number stiff;
	c_number r0;
	c_number rate;
	c_number r_ext;
	float3 centre;
};

void init_RepulsiveSphere_from_CPU(repulsive_sphere *cuda_force, RepulsiveSphere *cpu_force) {
	cuda_force->type = CUDA_REPULSIVE_SPHERE;
	cuda_force->stiff = cpu_force->_stiff;
	cuda_force->rate = cpu_force->_rate;
	cuda_force->r0 = cpu_force->_r0;
	cuda_force->r_ext = cpu_force->_r_ext;
	cuda_force->centre = make_float3(cpu_force->_center.x, cpu_force->_center.y, cpu_force->_center.z);
}

/**
 * @brief CUDA version of a RepulsiveSphereSmooth.
 */
struct repulsive_sphere_smooth {
	int type;
	c_number r0;
	c_number r_ext;
	c_number smooth;
	c_number alpha;
	c_number stiff;
	float3 centre;
};

void init_RepulsiveSphereSmooth_from_CPU(repulsive_sphere_smooth *cuda_force, RepulsiveSphereSmooth *cpu_force) {
	cuda_force->type = CUDA_REPULSIVE_SPHERE_SMOOTH;
	cuda_force->r0 = cpu_force->_r0;
	cuda_force->r_ext = cpu_force->_r_ext;
	cuda_force->smooth = cpu_force->_smooth;
	cuda_force->alpha = cpu_force->_alpha;
	cuda_force->stiff = cpu_force->_stiff;
	cuda_force->centre = make_float3(cpu_force->_center.x, cpu_force->_center.y, cpu_force->_center.z);
}

/**
 * @brief CUDA version of an LJWall.
 */
struct LJ_wall {
	int type;
	c_number stiff;
	c_number position;
	c_number cutoff;
	c_number sigma;
	int n;
	float3 dir;
};

void init_LJWall_from_CPU(LJ_wall *cuda_force, LJWall *cpu_force) {
	cuda_force->type = CUDA_LJ_WALL;
	cuda_force->stiff = cpu_force->_stiff;
	cuda_force->position = cpu_force->_position;
	cuda_force->n = cpu_force->_n;
	cuda_force->cutoff = cpu_force->_cutoff;
	cuda_force->sigma = cpu_force->_sigma;
	cuda_force->dir = make_float3(cpu_force->_direction.x, cpu_force->_direction.y, cpu_force->_direction.z);
}

/**
 * @brief CUDA version of a ConstantRateTorque.
 */
struct constant_rate_torque {
	int type;
	float3 center, pos0, axis, mask;
	c_number rate;
	c_number stiff;
	c_number F0;
};

void init_ConstantRateTorque_from_CPU(constant_rate_torque *cuda_force, ConstantRateTorque *cpu_force) {
	cuda_force->type = CUDA_CONSTANT_RATE_TORQUE;
	cuda_force->stiff = cpu_force->_stiff;
	cuda_force->F0 = cpu_force->_F0;
	cuda_force->rate = cpu_force->_rate;
	cuda_force->center = make_float3(cpu_force->_center.x, cpu_force->_center.y, cpu_force->_center.z);
	cuda_force->pos0 = make_float3(cpu_force->_pos0.x, cpu_force->_pos0.y, cpu_force->_pos0.z);
	cuda_force->axis = make_float3(cpu_force->_axis.x, cpu_force->_axis.y, cpu_force->_axis.z);
	cuda_force->mask = make_float3(cpu_force->_mask.x, cpu_force->_mask.y, cpu_force->_mask.z);
}

/**
 * @brief CUDA version of a GenericCentralForce.
 */
struct generic_constant_force {
	int type;
	c_number x, y, z;
	c_number F0;
	c_number inner_cut_off_sqr;
	c_number outer_cut_off_sqr;
};

void init_GenericCentralForce_from_CPU(generic_constant_force *cuda_force, GenericCentralForce *cpu_force) {
	cuda_force->type = CUDA_GENERIC_CENTRAL_FORCE;
	cuda_force->F0 = cpu_force->_F0;
	cuda_force->inner_cut_off_sqr = cpu_force->inner_cut_off_sqr;
	cuda_force->outer_cut_off_sqr = cpu_force->outer_cut_off_sqr;
	cuda_force->x = cpu_force->center.x;
	cuda_force->y = cpu_force->center.y;
	cuda_force->z = cpu_force->center.z;
}

/**
 * @brief CUDA version of an LJCone.
 */
struct LJ_cone {
	int type;
	float3 dir, pos0;
	c_number stiff;
	c_number cutoff;
	c_number sigma;
	c_number alpha;
	c_number sin_alpha;
	int n;
};

void init_LJCone_from_CPU(LJ_cone *cuda_force, LJCone *cpu_force) {
	cuda_force->type = CUDA_LJ_CONE;
	cuda_force->stiff = cpu_force->_stiff;
	cuda_force->n = cpu_force->_n;
	cuda_force->cutoff = cpu_force->_cutoff;
	cuda_force->sigma = cpu_force->_sigma;
	cuda_force->alpha = cpu_force->_alpha;
	cuda_force->sin_alpha = cpu_force->_sin_alpha;
	cuda_force->dir = make_float3(cpu_force->_direction.x, cpu_force->_direction.y, cpu_force->_direction.z);
	cuda_force->pos0 = make_float3(cpu_force->_pos0.x, cpu_force->_pos0.y, cpu_force->_pos0.z);
}

/**
 * @brief CUDA version of a RepulsiveEllipsoid.
 */
struct repulsive_ellipsoid {
	int type;
	c_number stiff;
	float3 r_1, r_2, centre;
};

void init_RepulsiveEllipsoid_from_CPU(repulsive_ellipsoid *cuda_force, RepulsiveEllipsoid *cpu_force) {
	cuda_force->type = CUDA_REPULSIVE_ELLIPSOID;
	cuda_force->stiff = cpu_force->_stiff;
	cuda_force->centre = make_float3(cpu_force->_centre.x, cpu_force->_centre.y, cpu_force->_centre.z);
	cuda_force->r_1 = make_float3(cpu_force->_r_1.x, cpu_force->_r_1.y, cpu_force->_r_1.z);
	cuda_force->r_2 = make_float3(cpu_force->_r_2.x, cpu_force->_r_2.y, cpu_force->_r_2.z);
}

/**
 * @brief CUDA version of a COMForce.
 */
struct COM_force {
    int type;
    c_number stiff;
    c_number r0;
    c_number rate;
    int n_com;
    int n_ref;
    int *com_indexes;
	int *ref_indexes;
};

void init_COMForce_from_CPU(COM_force *cuda_force, COMForce *cpu_force, bool first_time) {
	cuda_force->type = CUDA_COM_FORCE;
	cuda_force->stiff = cpu_force->_stiff;
	cuda_force->r0 = cpu_force->_r0;
	cuda_force->rate = cpu_force->_rate;
	cuda_force->n_com = cpu_force->_com_list.size();
	cuda_force->n_ref = cpu_force->_ref_list.size();

	std::vector<int> local_com_indexes;
	for(auto particle : cpu_force->_com_list) {
		local_com_indexes.push_back(particle->index);
	}

	std::vector<int> local_ref_indexes;
	for(auto particle : cpu_force->_ref_list){
		local_ref_indexes.push_back(particle->index);
	}

	if(first_time) {
		// TODO: this memory never gets free'd
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&cuda_force->com_indexes, sizeof(int) * local_com_indexes.size()));
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&cuda_force->ref_indexes, sizeof(int) * local_ref_indexes.size()));
	}

	CUDA_SAFE_CALL(cudaMemcpy(cuda_force->com_indexes, local_com_indexes.data(), sizeof(int) * local_com_indexes.size(), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(cuda_force->ref_indexes, local_ref_indexes.data(), sizeof(int) * local_ref_indexes.size(), cudaMemcpyHostToDevice));
}

/**
 * @brief CUDA version of a COMForce.
 */
struct lt_com_trap {
	int type;

	int p1a_size;
	int p2a_size;
    int *p1a;
    int *p2a;

    int xmin;
    int xmax;
    int N_grid;

    c_number dX;
    c_number *potential_grid;

    int mode;
	bool PBC;
};

void init_LTCOMTrap_from_CPU(lt_com_trap *cuda_force, LTCOMTrap *cpu_force, bool first_time) {
	cuda_force->type = CUDA_LR_COM_TRAP;
	cuda_force->xmin = cpu_force->xmin;
	cuda_force->xmax = cpu_force->xmax;
	cuda_force->N_grid = cpu_force->N_grid;
	cuda_force->dX = cpu_force->dX;
	cuda_force->mode = cpu_force->_mode;
	cuda_force->PBC = cpu_force->PBC;
	cuda_force->p1a_size = cpu_force->_p1a_ptr.size();
	cuda_force->p2a_size = cpu_force->_p2a_ptr.size();

	// copy the particle arrays
	std::vector<int> local_p1a;
	for(auto particle : cpu_force->_p1a_ptr) {
		local_p1a.push_back(particle->index);
	}
	std::vector<int> local_p2a;
	for(auto particle : cpu_force->_p2a_ptr) {
		local_p2a.push_back(particle->index);
	}

	// copy the potential array
	std::vector<c_number> local_grid;
	for(auto v : cpu_force->potential_grid) {
		local_grid.push_back(v);
	}

	if(first_time) {
		// TODO: this memory never gets free'd
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&cuda_force->p1a, sizeof(int) * local_p1a.size()));
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&cuda_force->p2a, sizeof(int) * local_p2a.size()));
		CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<c_number>(&cuda_force->potential_grid, sizeof(c_number) * local_grid.size()));
	}
	CUDA_SAFE_CALL(cudaMemcpy(cuda_force->p1a, local_p1a.data(), sizeof(int) * local_p1a.size(), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(cuda_force->p2a, local_p2a.data(), sizeof(int) * local_p2a.size(), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(cuda_force->potential_grid, local_grid.data(), sizeof(c_number) * local_grid.size(), cudaMemcpyHostToDevice));
}

/**
 * @brief Used internally by CUDA classes to provide an inheritance-like mechanism for external forces.
 */
union CUDA_trap {
	int type;
	constant_rate_force constant;
	mutual_trap mutual;
	moving_trap moving;
	lowdim_moving_trap lowdim;
	repulsion_plane repulsionplane;
	repulsion_plane_moving repulsionplanemoving;
	repulsive_sphere repulsivesphere;
	repulsive_sphere_smooth repulsivespheresmooth;
	LJ_wall ljwall;
	constant_rate_torque constantratetorque;
	generic_constant_force genericconstantforce;
	LJ_cone ljcone;
	repulsive_ellipsoid repulsiveellipsoid;
	COM_force comforce;
	lt_com_trap ltcomtrap;
};

#endif /* CUDAFORCES_H_ */
