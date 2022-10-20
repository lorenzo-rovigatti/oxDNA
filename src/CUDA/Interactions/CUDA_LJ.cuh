/* System constants */
__constant__ int MD_N[1];
__constant__ int MD_n_forces[1];
__constant__ float MD_sqr_rcut[3];
__constant__ float MD_sqr_sigma[3];
__constant__ float MD_epsilon[3];
__constant__ float MD_E_cut[3];
__constant__ int MD_LJ_n[3];

#include "../cuda_utils/CUDA_lr_common.cuh"

__device__ void _particle_particle_interaction(c_number4 &ppos, c_number4 &qpos, c_number4 &F, CUDABox *box) {
	int ptype = get_particle_type(ppos);
	int qtype = get_particle_type(qpos);
	int type = ptype + qtype;

	c_number4 r = box->minimum_image(ppos, qpos);

	c_number sqr_r = CUDA_DOT(r, r);

	c_number eps = MD_epsilon[type];
	c_number tmp = MD_sqr_sigma[type] / sqr_r;
	c_number lj_part = 1;
	for(int i = 0; i < MD_LJ_n[type] / 2; i++)
		lj_part *= tmp;
	c_number energy = 4.f * eps * (SQR(lj_part) - lj_part) - MD_E_cut[type];
	c_number force_module = -4.f * MD_LJ_n[type] * eps * (lj_part - 2.f * SQR(lj_part)) / sqr_r;

	if(sqr_r >= MD_sqr_rcut[type]) force_module = energy = (c_number) 0.f;

	F.x -= r.x * force_module;
	F.y -= r.y * force_module;
	F.z -= r.z * force_module;
	F.w += energy;
}

__global__ void lj_forces_edge(c_number4 *poss, c_number4 *forces, edge_bond *edge_list, int n_edges, CUDABox *box) {
	if(IND >= n_edges) return;

	c_number4 dF = make_c_number4(0.f, 0.f, 0.f, 0.f);

	edge_bond b = edge_list[IND];

	// get info for particle 1
	c_number4 ppos = poss[b.from];

	// get info for particle 2
	c_number4 qpos = poss[b.to];

	_particle_particle_interaction(ppos, qpos, dF, box);
	dF.w *= (c_number) 0.5f;

	int from_index = MD_N[0] * (IND % MD_n_forces[0]) + b.from;
	if((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z) > (c_number) 0.f) LR_atomicAddXYZ(forces + from_index, dF);

	dF.x = -dF.x;
	dF.y = -dF.y;
	dF.z = -dF.z;

	int to_index = MD_N[0] * (IND % MD_n_forces[0]) + b.to;
	if((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z) > (c_number) 0.f) LR_atomicAddXYZ(forces + to_index, dF);
}

__global__ void lj_forces(c_number4 *poss, c_number4 *forces, int *matrix_neighs, int *number_neighs, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 ppos = poss[IND];

	int num_neighs = NUMBER_NEIGHBOURS(IND, number_neighs);
	for(int j = 0; j < num_neighs; j++) {
		int k_index = NEXT_NEIGHBOUR(IND, j, matrix_neighs);

		if(k_index != IND) {
			c_number4 qpos = poss[k_index];
			_particle_particle_interaction(ppos, qpos, F, box);
		}
	}

	forces[IND] = F;
}
