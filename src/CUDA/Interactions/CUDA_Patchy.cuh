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

__device__ void _particle_particle_interaction(tmpnmbr &ppos, tmpnmbr &qpos, tmpnmbr &a1, tmpnmbr &a2, tmpnmbr &a3, tmpnmbr &b1, tmpnmbr &b2, tmpnmbr &b3, tmpnmbr &F, tmpnmbr &torque, CUDABox *box) {
	int ptype = get_particle_type(ppos);
	int qtype = get_particle_type(qpos);
	int type = ptype + qtype;

	tmpnmbr r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r);
	if(sqr_r >= MD_sqr_tot_rcut[type]) return;

	c_number part = powf(MD_sqr_sigma[type] / sqr_r, PATCHY_POWER * 0.5f);
	F.w += part - MD_E_cut[type];

	c_number force_module = PATCHY_POWER * part / sqr_r;
	F.x -= r.x * force_module;
	F.y -= r.y * force_module;
	F.z -= r.z * force_module;

	for(int pi = 0; pi < MD_N_patches[ptype]; pi++) {
		tmpnmbr ppatch = { a1.x * MD_base_patches[ptype][pi].x + a2.x * MD_base_patches[ptype][pi].y + a3.x * MD_base_patches[ptype][pi].z, a1.y * MD_base_patches[ptype][pi].x + a2.y * MD_base_patches[ptype][pi].y + a3.y * MD_base_patches[ptype][pi].z, a1.z * MD_base_patches[ptype][pi].x + a2.z * MD_base_patches[ptype][pi].y + a3.z * MD_base_patches[ptype][pi].z, 0 };
		ppatch *= MD_sigma[2 * ptype];
		for(int pj = 0; pj < MD_N_patches[qtype]; pj++) {
			tmpnmbr qpatch = { b1.x * MD_base_patches[qtype][pj].x + b2.x * MD_base_patches[qtype][pj].y + b3.x * MD_base_patches[qtype][pj].z, b1.y * MD_base_patches[qtype][pj].x + b2.y * MD_base_patches[qtype][pj].y + b3.y * MD_base_patches[qtype][pj].z, b1.z * MD_base_patches[qtype][pj].x + b2.z * MD_base_patches[qtype][pj].y + b3.z * MD_base_patches[qtype][pj].z, 0 };
			qpatch *= MD_sigma[2 * qtype];

			tmpnmbr patch_dist = { r.x + qpatch.x - ppatch.x, r.y + qpatch.y - ppatch.y, r.z + qpatch.z - ppatch.z, 0 };

			c_number dist = CUDA_DOT(patch_dist, patch_dist);
			if(dist < MD_sqr_patch_rcut[0]) {
				c_number r8b10 = dist * dist * dist * dist / MD_patch_pow_alpha[0];
				c_number exp_part = -1.001f * MD_epsilon[type] * expf(-(c_number) 0.5f * r8b10 * dist);

				F.w += exp_part - MD_patch_E_cut[type];

				tmpnmbr tmp_force = patch_dist * (5.f * exp_part * r8b10);

				torque -= _cross(ppatch, tmp_force);
				F.x -= tmp_force.x;
				F.y -= tmp_force.y;
				F.z -= tmp_force.z;
			}
		}
	}
}

// forces + second step without lists

__global__ void patchy_forces(tmpnmbr *poss, GPU_quat *orientations, tmpnmbr *forces, tmpnmbr *torques, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	tmpnmbr F = forces[IND];
	tmpnmbr T = make_tmpnmbr(0, 0, 0, 0);
	tmpnmbr ppos = poss[IND];
	GPU_quat po = orientations[IND];
	tmpnmbr a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND) {
			tmpnmbr qpos = poss[j];
			GPU_quat qo = orientations[j];
			get_vectors_from_quat(qo, b1, b2, b3);
			_particle_particle_interaction(ppos, qpos, a1, a2, a3, b1, b2, b3, F, T, box);
		}
	}

	T = _vectors_transpose_tmpnmbr_product(a1, a2, a3, T);

	forces[IND] = F;
	torques[IND] = T;
}

__global__ void patchy_forces_edge(tmpnmbr *poss, GPU_quat *orientations, tmpnmbr *forces, tmpnmbr *torques, edge_bond *edge_list, int n_edges, CUDABox *box) {
	if(IND >= n_edges) return;

	tmpnmbr dF = make_tmpnmbr(0, 0, 0, 0);
	tmpnmbr dT = make_tmpnmbr(0, 0, 0, 0);

	edge_bond b = edge_list[IND];

	// get info for particle 1
	tmpnmbr ppos = poss[b.from];
	GPU_quat po = orientations[b.from];

	// get info for particle 2
	tmpnmbr qpos = poss[b.to];
	GPU_quat qo = orientations[b.to];

	tmpnmbr a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);
	get_vectors_from_quat(qo, b1, b2, b3);

	_particle_particle_interaction(ppos, qpos, a1, a2, a3, b1, b2, b3, dF, dT, box);

	int from_index = MD_N[0] * (IND % MD_n_forces[0]) + b.from;
	if((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z) > (c_number) 0.f) LR_atomicAddXYZ(&(forces[from_index]), dF);
	if((dT.x * dT.x + dT.y * dT.y + dT.z * dT.z) > (c_number) 0.f) LR_atomicAddXYZ(&(torques[from_index]), _vectors_transpose_tmpnmbr_product(a1, a2, a3, dT));

	// Allen Eq. 6 pag 3:
	tmpnmbr dr = box->minimum_image(ppos, qpos); // returns qpos-ppos
	tmpnmbr crx = _cross(dr, dF);
	dT.x = -dT.x + crx.x;
	dT.y = -dT.y + crx.y;
	dT.z = -dT.z + crx.z;

	dF.x = -dF.x;
	dF.y = -dF.y;
	dF.z = -dF.z;

	int to_index = MD_N[0] * (IND % MD_n_forces[0]) + b.to;
	if((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z) > (c_number) 0.f) LR_atomicAddXYZ(&(forces[to_index]), dF);
	if((dT.x * dT.x + dT.y * dT.y + dT.z * dT.z) > (c_number) 0.f) LR_atomicAddXYZ(&(torques[to_index]), _vectors_transpose_tmpnmbr_product(b1, b2, b3, dT));
}
//Forces + second step with verlet lists

__global__ void patchy_forces(tmpnmbr *poss, GPU_quat *orientations, tmpnmbr *forces, tmpnmbr *torques, int *matrix_neighs, int *c_number_neighs, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	tmpnmbr F = forces[IND];
	tmpnmbr T = make_tmpnmbr(0, 0, 0, 0);
	tmpnmbr ppos = poss[IND];
	GPU_quat po = orientations[IND];
	tmpnmbr a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);

	int num_neighs = c_number_neighs[IND];

	for(int j = 0; j < num_neighs; j++) {
		int k_index = matrix_neighs[j * MD_N[0] + IND];

		tmpnmbr qpos = poss[k_index];
		GPU_quat qo = orientations[k_index];
		get_vectors_from_quat(qo, b1, b2, b3);
		_particle_particle_interaction(ppos, qpos, a1, a2, a3, b1, b2, b3, F, T, box);
	}

	T = _vectors_transpose_tmpnmbr_product(a1, a2, a3, T);

	forces[IND] = F;
	torques[IND] = T;
}
