/* System constants */
__constant__ int MD_N[1];
__constant__ int MD_n_forces[1];
__constant__ float MD_sqr_tot_rcut[3];
__constant__ float MD_sigma[3];
__constant__ float MD_sqr_sigma[3];
__constant__ float MD_epsilon[3];
__constant__ float MD_E_cut[3];
__constant__ float MD_sqr_patch_rcut[1];
__constant__ float MD_patch_pow_alpha[1];
__constant__ float MD_patch_E_cut[3];
__constant__ int MD_N_patches[2];
__constant__ float4 MD_base_patches[2][CUDA_MAX_PATCHES];
/* texture reference for positions */

#include "../cuda_utils/CUDA_lr_common.cuh"

__device__ void _particle_particle_interaction(c_number4 &ppos, c_number4 &qpos, c_number4 &a1, c_number4 &a2, c_number4 &a3, c_number4 &b1, c_number4 &b2, c_number4 &b3, c_number4 &F, c_number4 &torque, CUDABox *box) {
	int ptype = get_particle_type(ppos);
	int qtype = get_particle_type(qpos);
	int type = ptype + qtype;

	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r);
	if(sqr_r >= MD_sqr_tot_rcut[type]) return;

	c_number part = powf(MD_sqr_sigma[type] / sqr_r, PATCHY_POWER * 0.5f);
	F.w += part - MD_E_cut[type];

	c_number force_module = PATCHY_POWER * part / sqr_r;
	F.x -= r.x * force_module;
	F.y -= r.y * force_module;
	F.z -= r.z * force_module;

	for(int pi = 0; pi < MD_N_patches[ptype]; pi++) {
		c_number4 ppatch = { a1.x * MD_base_patches[ptype][pi].x + a2.x * MD_base_patches[ptype][pi].y + a3.x * MD_base_patches[ptype][pi].z, a1.y * MD_base_patches[ptype][pi].x + a2.y * MD_base_patches[ptype][pi].y + a3.y * MD_base_patches[ptype][pi].z, a1.z * MD_base_patches[ptype][pi].x + a2.z * MD_base_patches[ptype][pi].y + a3.z * MD_base_patches[ptype][pi].z, 0 };
		ppatch *= MD_sigma[2 * ptype];
		for(int pj = 0; pj < MD_N_patches[qtype]; pj++) {
			c_number4 qpatch = { b1.x * MD_base_patches[qtype][pj].x + b2.x * MD_base_patches[qtype][pj].y + b3.x * MD_base_patches[qtype][pj].z, b1.y * MD_base_patches[qtype][pj].x + b2.y * MD_base_patches[qtype][pj].y + b3.y * MD_base_patches[qtype][pj].z, b1.z * MD_base_patches[qtype][pj].x + b2.z * MD_base_patches[qtype][pj].y + b3.z * MD_base_patches[qtype][pj].z, 0 };
			qpatch *= MD_sigma[2 * qtype];

			c_number4 patch_dist = { r.x + qpatch.x - ppatch.x, r.y + qpatch.y - ppatch.y, r.z + qpatch.z - ppatch.z, 0 };

			c_number dist = CUDA_DOT(patch_dist, patch_dist);
			if(dist < MD_sqr_patch_rcut[0]) {
				c_number r8b10 = dist * dist * dist * dist / MD_patch_pow_alpha[0];
				c_number exp_part = -1.001f * MD_epsilon[type] * expf(-(c_number) 0.5f * r8b10 * dist);

				F.w += exp_part - MD_patch_E_cut[type];

				c_number4 tmp_force = patch_dist * (5.f * exp_part * r8b10);

				torque -= _cross(ppatch, tmp_force);
				F.x -= tmp_force.x;
				F.y -= tmp_force.y;
				F.z -= tmp_force.z;
			}
		}
	}
}

__global__ void patchy_forces_edge(c_number4 *poss, GPU_quat *orientations, c_number4 *forces, c_number4 *torques, edge_bond *edge_list, int n_edges, CUDABox *box) {
	if(IND >= n_edges) return;

	c_number4 dF = make_c_number4(0, 0, 0, 0);
	c_number4 dT = make_c_number4(0, 0, 0, 0);

	edge_bond b = edge_list[IND];

	// get info for particle 1
	c_number4 ppos = poss[b.from];
	GPU_quat po = orientations[b.from];

	// get info for particle 2
	c_number4 qpos = poss[b.to];
	GPU_quat qo = orientations[b.to];

	c_number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);
	get_vectors_from_quat(qo, b1, b2, b3);

	_particle_particle_interaction(ppos, qpos, a1, a2, a3, b1, b2, b3, dF, dT, box);

	int from_index = MD_N[0] * (IND % MD_n_forces[0]) + b.from;
	if((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z) > (c_number) 0.f) LR_atomicAddXYZ(&(forces[from_index]), dF);
	if((dT.x * dT.x + dT.y * dT.y + dT.z * dT.z) > (c_number) 0.f) LR_atomicAddXYZ(&(torques[from_index]), _vectors_transpose_c_number4_product(a1, a2, a3, dT));

	// Allen Eq. 6 pag 3:
	c_number4 dr = box->minimum_image(ppos, qpos); // returns qpos-ppos
	c_number4 crx = _cross(dr, dF);
	dT.x = -dT.x + crx.x;
	dT.y = -dT.y + crx.y;
	dT.z = -dT.z + crx.z;

	dF.x = -dF.x;
	dF.y = -dF.y;
	dF.z = -dF.z;

	int to_index = MD_N[0] * (IND % MD_n_forces[0]) + b.to;
	if((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z) > (c_number) 0.f) LR_atomicAddXYZ(&(forces[to_index]), dF);
	if((dT.x * dT.x + dT.y * dT.y + dT.z * dT.z) > (c_number) 0.f) LR_atomicAddXYZ(&(torques[to_index]), _vectors_transpose_c_number4_product(b1, b2, b3, dT));
}

__global__ void patchy_forces(c_number4 *poss, GPU_quat *orientations, c_number4 *forces, c_number4 *torques, int *matrix_neighs, int *number_neighs, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = make_c_number4(0, 0, 0, 0);
	c_number4 ppos = poss[IND];
	GPU_quat po = orientations[IND];
	c_number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);

	int num_neighs = NUMBER_NEIGHBOURS(IND, number_neighs);
	for(int j = 0; j < num_neighs; j++) {
		int k_index = NEXT_NEIGHBOUR(IND, j, matrix_neighs);

		if(k_index != IND) {
			c_number4 qpos = poss[k_index];
			GPU_quat qo = orientations[k_index];
			get_vectors_from_quat(qo, b1, b2, b3);
			_particle_particle_interaction(ppos, qpos, a1, a2, a3, b1, b2, b3, F, T, box);
		}
	}

	T = _vectors_transpose_c_number4_product(a1, a2, a3, T);

	forces[IND] = F;
	torques[IND] = T;
}
