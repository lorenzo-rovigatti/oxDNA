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
#define CUDA_REPULSIVE_ELLIPSOID 12

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
};

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

/**
 * @brief CUDA version of a RepulsionPlane.
 */

struct repulsion_plane {
	int type;
	c_number stiff;
	c_number position;
	float3 dir;
};

/**
 * @brief CUDA version of a RepulsionPlaneMoving.
 */

struct repulsion_plane_moving {
	int type;
	c_number stiff;
	int low_idx, high_idx;
	float3 dir;
};

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

/**
 * @brief CUDA version of a RepulsiveEllipsoid.
 */
struct repulsive_ellipsoid {
	int type;
	c_number stiff;
	float3 r_1, r_2, centre;
};

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
};

#endif /* CUDAFORCES_H_ */
