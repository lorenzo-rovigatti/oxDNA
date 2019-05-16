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
#define CUDA_REPULSIVE_SPHERE 7
#define CUDA_CONSTANT_RATE_TORQUE 8
#define CUDA_GENERIC_CENTRAL_FORCE 9
#define CUDA_LJ_CONE 10
#define CUDA_REPULSIVE_SPHERE_SMOOTH 11

/**
 * @brief CUDA version of a ConstantRateForce.
 */
template<typename number>
struct constant_rate_force {
	int type;
	number x, y, z;
	bool dir_as_centre;
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
	number rate;
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
	int low_idx, high_idx;
	float3 dir;
};

/**
 * @brief CUDA version of a RepulsiveSphere.
 */
template<typename number>
struct repulsive_sphere {
	int type;
	number stiff;
	number r0;
	number rate;
	number r_ext;
	float3 centre;
};

/**
 * @brief CUDA version of a RepulsiveSphereSmooth.
 */
template<typename number>
struct repulsive_sphere_smooth {
	int type;
	number r0;
	number r_ext;
	number smooth;
	number alpha;
	number stiff;
	float3 centre;
};

/**
 * @brief CUDA version of an LJWall.
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
 * @brief CUDA version of a ConstantRateTorque.
 */
template<typename number>
struct constant_rate_torque {
	int type;
	float3 center, pos0, axis, mask;
	number rate;
	number stiff;
	number F0;
};

/**
 * @brief CUDA version of a GenericCentralForce.
 */
template<typename number>
struct generic_constant_force {
	int type;
	number x, y, z;
	number F0;
	number inner_cut_off_sqr;
	number outer_cut_off_sqr;
};

/**
 * @brief CUDA version of an LJCone.
 */
template<typename number>
struct LJ_cone {
	int type;
	float3 dir, pos0;
	number stiff;
	number cutoff;
	number sigma;
	number alpha;
	number sin_alpha;
	int n;
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
	repulsive_sphere<number> repulsivesphere;
	repulsive_sphere_smooth<number> repulsivespheresmooth;
	LJ_wall<number> ljwall;
	constant_rate_torque<number> constantratetorque;
	generic_constant_force<number> genericconstantforce;
	LJ_cone<number> ljcone;
};

#endif /* CUDAFORCES_H_ */
