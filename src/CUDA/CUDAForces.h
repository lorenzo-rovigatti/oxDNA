/*
 * CUDAForces.h
 *
 *  Created on: 26/giu/2013
 *      Author: lorenzo
 */

#ifndef CUDAFORCES_H_
#define CUDAFORCES_H_

#define CUDA_TRAP_NO_FORCE -1
#define CUDA_TRAP_CONSTANT 0
#define CUDA_TRAP_MUTUAL 1
#define CUDA_TRAP_MOVING 2
#define CUDA_REPULSION_PLANE 3
#define CUDA_REPULSION_PLANE_MOVING 4
#define CUDA_TRAP_MOVING_LOWDIM 5
#define CUDA_LJ_WALL 6

/**
 * @brief CUDA version of a ConstantRateForce.
 */
template<typename number>
struct constant_rate_force {
	int type;
	number x, y, z;
	number rate;
	number F0;
};

/**
 * @brief CUDA version of a MutualTrap.
 */
template<typename number>
struct mutual_trap {
	int type;
	number stiff;
	number r0;
	int p_ind;
	bool PBC;
};

/**
 * @brief CUDA version of a MovingTrap.
 */
template<typename number>
struct moving_trap {
	int type;
	number stiff;
	number rate;
	float3 pos0;
	float3 dir;
};

/**
 * @brief CUDA version of a LowdimMovingTrap.
 */
template<typename number>
struct lowdim_moving_trap {
	int type;
	number stiff;
	number rate;
	float3 pos0;
	float3 dir;
	bool visX;
	bool visY;
	bool visZ;
};

/**
 * @brief CUDA version of a RepulsionPlane.
 */
template<typename number>
struct repulsion_plane {
	int type;
	number stiff;
	number position;
	float3 dir;
};

/**
 * @brief CUDA version of a RepulsionPlaneMoving.
 */
template<typename number>
struct repulsion_plane_moving {
	int type;
	number stiff;
	int p_ind;
	float3 dir;
};

/**
 * @brief CUDA version of a LJWall.
 */
template<typename number>
struct LJ_wall {
	int type;
	number stiff;
	number position;
	number cutoff;
	number sigma;
	int n;
	float3 dir;
};

/**
 * @brief Used internally by CUDA classes to provide an inheritance-like mechanism for external forces.
 */
template<typename number>
union CUDA_trap {
	int type;
	constant_rate_force<number> constant;
	mutual_trap<number> mutual;
	moving_trap<number> moving;
	lowdim_moving_trap<number> lowdim;
	repulsion_plane<number> repulsionplane;
	repulsion_plane_moving<number> repulsionplanemoving;
	LJ_wall<number> ljwall;
};

#endif /* CUDAFORCES_H_ */
