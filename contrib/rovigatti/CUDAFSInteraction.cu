/*
 * CUDAFSInteraction.cu
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include "CUDAFSInteraction.h"

#include "../Lists/CUDASimpleVerletList.h"
#include "../Lists/CUDANoList.h"

/* BEGIN CUDA */

/* System constants */
__constant__ int MD_N[1];
__constant__ int MD_n_forces[1];
__constant__ int MD_N_patches[2];
__constant__ bool MD_one_component[1];
__constant__ float MD_box_side[1];
__constant__ float MD_sqr_rcut[1];
__constant__ float MD_sqr_rep_rcut[1];
__constant__ float MD_sqr_patch_rcut[1];
__constant__ float MD_sigma_ss[1];
__constant__ float MD_rcut_ss[1];
__constant__ float MD_lambda[1];
__constant__ float MD_A_part[1], MD_B_part[1];
__constant__ float4 MD_base_patches[2][CUDA_MAX_FS_PATCHES];

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
__device__ void _particle_particle_interaction(number4 &ppos, number4 &qpos, number4 &a1, number4 &a2, number4 &a3, number4 &b1, number4 &b2, number4 &b3, number4 &F, number4 &torque) {
	int ptype = get_particle_type<number, number4>(ppos);
	int qtype = get_particle_type<number, number4>(qpos);

	number4 r = minimum_image<number, number4>(ppos, qpos);
	number sqr_r = CUDA_DOT(r, r);
	if(sqr_r >= MD_sqr_rcut[0]) return;

	// centre-centre
	number ir2 = 1.f / sqr_r;
	number lj_part = ir2*ir2*ir2;
	number force_module = -24.f * (lj_part - 2.f*SQR(lj_part)) / sqr_r;
	if(sqr_r >= MD_sqr_rep_rcut[0]) force_module = 0.f;
	F.x -= r.x * force_module;
	F.y -= r.y * force_module;
	F.z -= r.z * force_module;

	// TODO may be improved by removing branching
	if(ptype == qtype && !MD_one_component[0]) return;

	for(int pi = 0; pi < MD_N_patches[ptype]; pi++) {
		number4 ppatch = {
			a1.x*MD_base_patches[ptype][pi].x + a2.x*MD_base_patches[ptype][pi].y + a3.x*MD_base_patches[ptype][pi].z,
			a1.y*MD_base_patches[ptype][pi].x + a2.y*MD_base_patches[ptype][pi].y + a3.y*MD_base_patches[ptype][pi].z,
			a1.z*MD_base_patches[ptype][pi].x + a2.z*MD_base_patches[ptype][pi].y + a3.z*MD_base_patches[ptype][pi].z,
			0
		};

		for(int pj = 0; pj < MD_N_patches[qtype]; pj++) {
			number4 qpatch = {
				b1.x*MD_base_patches[qtype][pj].x + b2.x*MD_base_patches[qtype][pj].y + b3.x*MD_base_patches[qtype][pj].z,
				b1.y*MD_base_patches[qtype][pj].x + b2.y*MD_base_patches[qtype][pj].y + b3.y*MD_base_patches[qtype][pj].z,
				b1.z*MD_base_patches[qtype][pj].x + b2.z*MD_base_patches[qtype][pj].y + b3.z*MD_base_patches[qtype][pj].z,
				0
			};

			number4 patch_dist = {
				r.x + qpatch.x - ppatch.x,
				r.y + qpatch.y - ppatch.y,
				r.z + qpatch.z - ppatch.z,
				0
			};

			number dist = CUDA_DOT(patch_dist, patch_dist);
			if(dist < MD_sqr_patch_rcut[0]) {
				number r_p = sqrtf(dist);
				number exp_part = expf(MD_sigma_ss[0] / (r_p - MD_rcut_ss[0]));
				number energy_part = MD_A_part[0] * exp_part * (MD_B_part[0]/SQR(dist) - 1.);

				number force_mod  = MD_A_part[0] * exp_part * (4.*MD_B_part[0]/(SQR(dist)*r_p)) + MD_sigma_ss[0] * energy_part / SQR(r_p - MD_rcut_ss[0]);
				number4 tmp_force = patch_dist * (force_mod / r_p);

				torque -= _cross<number, number4>(ppatch, tmp_force);
				F.x -= tmp_force.x;
				F.y -= tmp_force.y;
				F.z -= tmp_force.z;
			}
		}
	}
}

// forces + second step without lists
template <typename number, typename number4>
__global__ void FS_forces(number4 *poss, GPU_quat<number> *orientations, number4 *forces, number4 *torques) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	number4 T = make_number4<number, number4>(0, 0, 0, 0);
	number4 ppos = poss[IND];
	GPU_quat<number> po = orientations[IND];
	number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat<number,number4>(po, a1, a2, a3);

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND) {
			number4 qpos = poss[j];
			GPU_quat<number> qo = orientations[j];
			get_vectors_from_quat<number,number4>(qo, b1, b2, b3);
			_particle_particle_interaction<number, number4>(ppos, qpos, a1, a2, a3, b1, b2, b3, F, T);
		}
	}

	T = _vectors_transpose_number4_product(a1, a2, a3, T);

	forces[IND] = F;
	torques[IND] = T;
}

template <typename number, typename number4>
__global__ void FS_forces_edge(number4 *poss, GPU_quat<number> *orientations, number4 *forces, number4 *torques, edge_bond *edge_list,  int n_edges) {
	if(IND >= n_edges) return;

	number4 dF = make_number4<number, number4>(0, 0, 0, 0);
	number4 dT = make_number4<number, number4>(0, 0, 0, 0);

	edge_bond b = edge_list[IND];

	// get info for particle 1
	number4 ppos = poss[b.from];
	GPU_quat<number> po = orientations[b.from];

	// get info for particle 2
	number4 qpos = poss[b.to];
	GPU_quat<number> qo = orientations[b.to];

	number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat<number,number4>(po, a1, a2, a3);
	get_vectors_from_quat<number,number4>(qo, b1, b2, b3);

	_particle_particle_interaction<number, number4>(ppos, qpos, a1, a2, a3, b1, b2, b3, dF, dT);

	int from_index = MD_N[0]*(IND % MD_n_forces[0]) + b.from;
	if((dF.x*dF.x + dF.y*dF.y + dF.z*dF.z) > (number)0.f) LR_atomicAddXYZ(&(forces[from_index]), dF);
	if((dT.x*dT.x + dT.y*dT.y + dT.z*dT.z) > (number)0.f) LR_atomicAddXYZ(&(torques[from_index]), _vectors_transpose_number4_product(a1, a2, a3, dT));

	// Allen Eq. 6 pag 3:
	number4 dr = minimum_image<number, number4>(ppos, qpos); // returns qpos-ppos
	number4 crx = _cross<number, number4>(dr, dF);
	dT.x = -dT.x + crx.x;
	dT.y = -dT.y + crx.y;
	dT.z = -dT.z + crx.z;

	dF.x = -dF.x;
	dF.y = -dF.y;
	dF.z = -dF.z;

	int to_index = MD_N[0]*(IND % MD_n_forces[0]) + b.to;
	if((dF.x*dF.x + dF.y*dF.y + dF.z*dF.z) > (number)0.f) LR_atomicAddXYZ(&(forces[to_index]), dF);
	if((dT.x*dT.x + dT.y*dT.y + dT.z*dT.z) > (number)0.f) LR_atomicAddXYZ(&(torques[to_index]), _vectors_transpose_number4_product(b1, b2, b3, dT));
}

//Forces + second step with verlet lists
template <typename number, typename number4>
__global__ void FS_forces(number4 *poss, GPU_quat<number> *orientations, number4 *forces, number4 *torques, int *matrix_neighs, int *number_neighs) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	number4 T = make_number4<number, number4>(0, 0, 0, 0);
	number4 ppos = poss[IND];
	GPU_quat<number> po = orientations[IND];
	number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat<number,number4>(po, a1, a2, a3);

	int num_neighs = number_neighs[IND];

	for(int j = 0; j < num_neighs; j++) {
		int k_index = matrix_neighs[j*MD_N[0] + IND];

		number4 qpos = poss[k_index];
		GPU_quat<number> qo = orientations[k_index];
		get_vectors_from_quat<number,number4>(qo, b1, b2, b3);
		_particle_particle_interaction<number, number4>(ppos, qpos, a1, a2, a3, b1, b2, b3, F, T);
	}

	T = _vectors_transpose_number4_product(a1, a2, a3, T);

	forces[IND] = F;
	torques[IND] = T;
}

/* END CUDA PART */

#define HALF_ISQRT3 0.28867513459481292f

template<typename number, typename number4>
CUDAFSInteraction<number, number4>::CUDAFSInteraction() : CUDABaseInteraction<number, number4>(), FSInteraction<number>() {

}

template<typename number, typename number4>
CUDAFSInteraction<number, number4>::~CUDAFSInteraction() {

}

template<typename number, typename number4>
void CUDAFSInteraction<number, number4>::get_settings(input_file &inp) {
	FSInteraction<number>::get_settings(inp);
}

template<typename number, typename number4>
void CUDAFSInteraction<number, number4>::cuda_init(number box_side, int N) {
	CUDABaseInteraction<number, number4>::cuda_init(box_side, N);
	FSInteraction<number>::init();

	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &N, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_one_component, &this->_one_component, sizeof(bool)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N_patches, &this->_N_patches, sizeof(int)) );
	if(!this->_one_component) CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N_patches, &this->_N_patches_B, sizeof(int), sizeof(int)) );

	float f_copy = box_side;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_box_side, &f_copy, sizeof(float)) );
	f_copy = this->_sqr_rcut;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_sqr_rcut, &f_copy, sizeof(float)) );
	f_copy = this->_sqr_rep_rcut;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_sqr_rep_rcut, &f_copy, sizeof(float)) );
	f_copy = this->_sqr_patch_rcut;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_sqr_patch_rcut, &f_copy, sizeof(float)) );
	f_copy = this->_sigma_ss;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_sigma_ss, &f_copy, sizeof(float)) );
	f_copy = this->_rcut_ss;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_rcut_ss, &f_copy, sizeof(float)) );
	f_copy = this->_lambda;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_lambda, &f_copy, sizeof(float)) );
	f_copy = this->_A_part;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_A_part, &f_copy, sizeof(float)) );
	f_copy = this->_B_part;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_B_part, &f_copy, sizeof(float)) );

	float4 base_patches[CUDA_MAX_FS_PATCHES];

	// ugly...
	int limit = (this->_one_component) ? 1 : 2;
	int n_patches = this->_N_patches;
	for(int i = 0; i < limit; i++) {
		switch(n_patches) {
		case 2: {
			base_patches[0] = make_float4(0, 0.5, 0, 0);
			base_patches[1] = make_float4(0, -0.5, 0, 0);
			break;
		}
		case 3: {
			number cos120 = cos(2 * M_PI / 3.);
			number sin120 = sin(2 * M_PI / 3.);

			base_patches[0] = make_float4(0, 1, 0, 0);
			base_patches[1] = make_float4(cos120, -sin120, 0, 0);
			base_patches[2] = make_float4(-cos120, -sin120, 0, 0);
			break;
		}
		case 4: {
			base_patches[0] = make_float4(-HALF_ISQRT3, -HALF_ISQRT3,  HALF_ISQRT3, 0);
			base_patches[1] = make_float4( HALF_ISQRT3, -HALF_ISQRT3, -HALF_ISQRT3, 0);
			base_patches[2] = make_float4( HALF_ISQRT3,  HALF_ISQRT3,  HALF_ISQRT3, 0);
			base_patches[3] = make_float4(-HALF_ISQRT3,  HALF_ISQRT3, -HALF_ISQRT3, 0);
			break;
		}
		default:
			throw oxDNAException("Unsupported number of patches %d", n_patches);
		}

		for(int j = 0; j < n_patches; j++) {
			number factor = 0.5 / sqrt(CUDA_DOT(base_patches[j], base_patches[j]));
			base_patches[j].x *= factor;
			base_patches[j].y *= factor;
			base_patches[j].z *= factor;
		}

		// fourth argument is the offset
		CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_base_patches, base_patches, sizeof(float4)*n_patches, i*sizeof(float4)*CUDA_MAX_FS_PATCHES) );
		n_patches = this->_N_patches_B;
	}

	if(this->_N_patches > CUDA_MAX_FS_PATCHES) throw oxDNAException("CUDA supports only particles with up to %d patches", CUDA_MAX_FS_PATCHES);
	if(this->_use_edge) CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_n_forces, &this->_n_forces, sizeof(int)) );
}

template<typename number, typename number4>
void CUDAFSInteraction<number, number4>::compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds) {
	CUDASimpleVerletList<number, number4> *_v_lists = dynamic_cast<CUDASimpleVerletList<number, number4> *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) {
				FS_forces_edge<number, number4>
					<<<(_v_lists->_N_edges - 1)/(this->_launch_cfg.threads_per_block) + 1, this->_launch_cfg.threads_per_block>>>
					//(d_poss, d_orientations, d_forces, d_torques, _v_lists->_d_edge_list, _v_lists->_N_edges);
					(d_poss, d_orientations, this->_d_edge_forces, this->_d_edge_torques, _v_lists->_d_edge_list, _v_lists->_N_edges);
				CUT_CHECK_ERROR("forces_second_step FS forces_edge");

				this->_sum_edge_forces_torques(d_forces, d_torques);
			}
			else {
				FS_forces<number, number4>
					<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
					(d_poss, d_orientations, d_forces, d_torques, _v_lists->_d_matrix_neighs, _v_lists->_d_number_neighs);
				CUT_CHECK_ERROR("forces_second_step FS simple_lists error");
			}
	}

	CUDANoList<number, number4> *_no_lists = dynamic_cast<CUDANoList<number, number4> *>(lists);
	if(_no_lists != NULL) {
		FS_forces<number, number4>
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_orientations, d_forces, d_torques);
		CUT_CHECK_ERROR("forces_second_step FS no_lists error");
	}
}

template class CUDAFSInteraction<float, float4>;
template class CUDAFSInteraction<double, LR_double4>;
