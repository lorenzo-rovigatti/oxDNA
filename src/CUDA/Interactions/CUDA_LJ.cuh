/* System constants */
__constant__ int MD_N[1];
__constant__ int MD_n_forces[1];
__constant__ float MD_box_side[1];
__constant__ float MD_sqr_rcut[3];
__constant__ float MD_sqr_sigma[3];
__constant__ float MD_epsilon[3];
__constant__ int MD_LJ_n[3];

#include "../cuda_utils/CUDA_lr_common.cuh"

template <typename number, typename number4>
__device__ number4 minimum_image(number4 &r_i, number4 &r_j) {
	number dx = r_j.x - r_i.x;
	number dy = r_j.y - r_i.y;
	number dz = r_j.z - r_i.z;

	dx -= floorf(dx/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];
	dy -= floorf(dy/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];
	dz -= floorf(dz/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];

	return make_number4<number, number4>(dx, dy, dz, (number) 0.f);
}

template <typename number, typename number4>
__device__ void _particle_particle_interaction(number4 &ppos, number4 &qpos, number4 &F) {
	int ptype = get_particle_type<number, number4>(ppos);
	int qtype = get_particle_type<number, number4>(qpos);
	int type = ptype + qtype;

	number4 r = minimum_image<number, number4>(ppos, qpos);

	number sqr_r = CUDA_DOT(r, r);

	number eps = MD_epsilon[type];
	number tmp = MD_sqr_sigma[type] / sqr_r;
	number lj_part = 1;
	for(int i = 0; i < MD_LJ_n[type]/2; i++) lj_part *= tmp;
	number force_module = -4.f * MD_LJ_n[type] * eps * (lj_part - 2.f*SQR(lj_part)) / sqr_r;

	if(sqr_r >= MD_sqr_rcut[type]) force_module = (number) 0.f;

	F.x -= r.x * force_module;
	F.y -= r.y * force_module;
	F.z -= r.z * force_module;
}

// forces + second step without lists
template <typename number, typename number4>
__global__ void lj_forces(number4 *poss, number4 *forces) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	number4 ppos = poss[IND];

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND) {
			number4 qpos = poss[j];

			_particle_particle_interaction<number, number4>(ppos, qpos, F);
		}
	}

	forces[IND] = F;
}

template <typename number, typename number4>
__global__ void lj_forces_edge(number4 *poss, number4 *forces, edge_bond *edge_list, int n_edges) {
	if(IND >= n_edges) return;

	number4 dF = make_number4<number, number4>(0.f, 0.f, 0.f, 0.f);

	edge_bond b = edge_list[IND];

	// get info for particle 1
	number4 ppos = poss[b.from];

	// get info for particle 2
	number4 qpos = poss[b.to];

	_particle_particle_interaction<number, number4>(ppos, qpos, dF);
	dF.w *= (number) 0.5f;

	int from_index = MD_N[0]*(IND % MD_n_forces[0]) + b.from;
	if((dF.x*dF.x + dF.y*dF.y + dF.z*dF.z) > (number)0.f) LR_atomicAddXYZ(forces + from_index, dF);

	dF.x = -dF.x;
	dF.y = -dF.y;
	dF.z = -dF.z;

	int to_index = MD_N[0]*(IND % MD_n_forces[0]) + b.to;
	if ((dF.x*dF.x + dF.y*dF.y + dF.z*dF.z) > (number)0.f) LR_atomicAddXYZ(forces + to_index, dF);
}

// forces + second step with verlet lists
template <typename number, typename number4>
__global__ void lj_forces(number4 *poss, number4 *forces, int *matrix_neighs, int *number_neighs) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	number4 ppos = poss[IND];

	int num_neighs = number_neighs[IND];

	for(int j = 0; j < num_neighs; j++) {
		int k_index = matrix_neighs[j*MD_N[0] + IND];

		number4 qpos = poss[k_index];
		_particle_particle_interaction<number, number4>(ppos, qpos, F);
	}


	forces[IND] = F;
}
