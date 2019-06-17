/*
 * CUDAPolymerInteraction.cu
 *
 *  Created on: 08/apr/2019
 *      Author: lorenzo
 */

#include "CUDAPolymerInteraction.h"

#include "CUDA/Lists/CUDASimpleVerletList.h"
#include "CUDA/Lists/CUDANoList.h"

/* BEGIN CUDA */

/* System constants */
__constant__ int MD_N[1];
__constant__ int MD_n_forces[1];

__constant__ float MD_sqr_rfene[1];
__constant__ float MD_sqr_rcut[1];
__constant__ float MD_sqr_rep_rcut[1];
__constant__ float MD_Polymer_lambda[1];

#include "CUDA/cuda_utils/CUDA_lr_common.cuh"

template <typename number, typename number4>
__device__ void _nonbonded(number4 &r, int int_type, number4 &F) {
	number sqr_r = CUDA_DOT(r, r);

	number energy = 0.f;
	// this number is the module of the force over r, so we don't have to divide the distance
	// vector for its module
	number force_mod = 0.f;

	if(sqr_r < MD_sqr_rep_rcut[0]) {
		number part = powf(1.f / sqr_r, 3.f);
		energy += 4.f * part * (part - 1.f) + 1.f;
		force_mod += 24.f * part * (2.f*part - 1.f) / sqr_r;
	}

	if(int_type == 2) {
		if(sqr_r < MD_sqr_rep_rcut[0]) energy -= MD_Polymer_lambda[0];
		else {
			number part = powf(1.f/sqr_r, 3.f);
			energy += 24.f * part * (part - 1.f);
			force_mod += 24.f * MD_Polymer_lambda[0] * part*(2*part - 1.f) / sqr_r;
		}
	}

	if(sqr_r > MD_sqr_rcut[0]) energy = force_mod = (number) 0.f;

	F.x -= r.x*force_mod;
	F.y -= r.y*force_mod;
	F.z -= r.z*force_mod;
	F.w += energy;
}

template<typename number, typename number4>
__device__ void _fene(number4 &r, number4 &F) {
	number sqr_r = CUDA_DOT(r, r);

	number energy = -15.f * MD_sqr_rfene[0] * logf(1.f - sqr_r / MD_sqr_rfene[0]);

	// this number is the module of the force over r, so we don't have to divide the distance
	// vector by its module
	number force_mod = -30.f * MD_sqr_rfene[0] / (MD_sqr_rfene[0] - sqr_r);
	F.x -= r.x * force_mod;
	F.y -= r.y * force_mod;
	F.z -= r.z * force_mod;
	F.w += energy;
}

template <typename number, typename number4>
__device__ void _particle_particle_bonded_interaction(number4 &ppos, number4 &qpos, number4 &F) {
	number4 r = qpos - ppos;
	// bonded interactions are purely repulsive, so we set the int_type to 0 for all pairs
	_nonbonded<number, number4>(r, 0, F);
	_fene<number, number4>(r, F);
}

template <typename number, typename number4>
__device__ void _particle_particle_interaction(number4 &ppos, number4 &qpos, number4 &F, CUDABox<number, number4> *box) {
	int ptype = get_particle_type<number, number4>(ppos);
	int qtype = get_particle_type<number, number4>(qpos);
	int int_type = ptype + qtype;

	number4 r = box->minimum_image(ppos, qpos);
	_nonbonded<number, number4>(r, int_type, F);
}

// forces + second step without lists
template <typename number, typename number4>
__global__ void polymer_forces(number4 *poss, number4 *forces, LR_bonds *bonds, CUDABox<number, number4> *box) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	LR_bonds bs = bonds[IND];
	number4 ppos = poss[IND];

	if(bs.n3 != P_INVALID) {
		number4 qpos = poss[bs.n3];
		_particle_particle_bonded_interaction<number, number4>(ppos, qpos, F);
	}

	if(bs.n5 != P_INVALID) {
		number4 qpos = poss[bs.n5];
		_particle_particle_bonded_interaction<number, number4>(ppos, qpos, F);
	}

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND && bs.n3 != j && bs.n5 != j) {
			number4 qpos = poss[j];
			_particle_particle_interaction<number, number4>(ppos, qpos, F, box);
		}
	}

	forces[IND] = F;
}

template <typename number, typename number4>
__global__ void polymer_forces_edge_nonbonded(number4 *poss, number4 *forces, edge_bond *edge_list, int n_edges, CUDABox<number, number4> *box) {
	if(IND >= n_edges) return;

	number4 dF = make_number4<number, number4>(0, 0, 0, 0);

	edge_bond b = edge_list[IND];

	// get info for particle 1
	number4 ppos = poss[b.from];

	// get info for particle 2
	number4 qpos = poss[b.to];

	_particle_particle_interaction<number, number4>(ppos, qpos, dF, box);

	dF.w *= (number) 0.5f;

	int from_index = MD_N[0]*(IND % MD_n_forces[0]) + b.from;
	//int from_index = MD_N[0]*(b.n_from % MD_n_forces[0]) + b.from;
	if((dF.x*dF.x + dF.y*dF.y + dF.z*dF.z + dF.w*dF.w) > (number)0.f) LR_atomicAddXYZ(&(forces[from_index]), dF);

	// Allen Eq. 6 pag 3:
	number4 dr = box->minimum_image(ppos, qpos); // returns qpos-ppos
	number4 crx = _cross<number, number4> (dr, dF);

	dF.x = -dF.x;
	dF.y = -dF.y;
	dF.z = -dF.z;

	int to_index = MD_N[0]*(IND % MD_n_forces[0]) + b.to;
	if((dF.x*dF.x + dF.y*dF.y + dF.z*dF.z + dF.w*dF.w) > (number)0.f) LR_atomicAddXYZ(&(forces[to_index]), dF);
}

// bonded interactions for edge-based approach
template <typename number, typename number4>
__global__ void polymer_forces_edge_bonded(number4 *poss, number4 *forces, LR_bonds *bonds) {
	if(IND >= MD_N[0]) return;

	number4 F0;

	F0.x = forces[IND].x;
	F0.y = forces[IND].y;
	F0.z = forces[IND].z;
	F0.w = forces[IND].w;

	number4 dF = make_number4<number, number4>(0, 0, 0, 0);
	number4 ppos = poss[IND];
	LR_bonds bs = bonds[IND];

	if(bs.n3 != P_INVALID) {
		number4 qpos = poss[bs.n3];
		_particle_particle_bonded_interaction<number, number4>(ppos, qpos, dF);
	}
	if(bs.n5 != P_INVALID) {
		number4 qpos = poss[bs.n5];
		_particle_particle_bonded_interaction<number, number4>(ppos, qpos, dF);
	}

	forces[IND] = (dF + F0);
}

// forces + second step with verlet lists
template <typename number, typename number4>
__global__ void polymer_forces(number4 *poss, number4 *forces, int *matrix_neighs, int *number_neighs, LR_bonds *bonds, CUDABox<number, number4> *box) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	number4 ppos = poss[IND];
	LR_bonds bs = bonds[IND];

	if(bs.n3 != P_INVALID) {
		number4 qpos = poss[bs.n3];
		_particle_particle_bonded_interaction<number, number4>(ppos, qpos, F);
	}
	if(bs.n5 != P_INVALID) {
		number4 qpos = poss[bs.n5];
		_particle_particle_bonded_interaction<number, number4>(ppos, qpos, F);
	}

	const int num_neighs = number_neighs[IND];
	for(int j = 0; j < num_neighs; j++) {
		const int k_index = matrix_neighs[j*MD_N[0] + IND];

		number4 qpos = poss[k_index];
		_particle_particle_interaction<number, number4>(ppos, qpos, F, box);
	}

	forces[IND] = F;
}

/* END CUDA */

template<typename number, typename number4>
CUDAPolymerInteraction<number, number4>::CUDAPolymerInteraction() {

}

template<typename number, typename number4>
CUDAPolymerInteraction<number, number4>::~CUDAPolymerInteraction() {

}

template<typename number, typename number4>
void CUDAPolymerInteraction<number, number4>::get_settings(input_file &inp) {
	PolymerInteraction<number>::get_settings(inp);

	int sort_every;
	if(getInputInt(&inp, "CUDA_sort_every", &sort_every, 0) == KEY_FOUND) {
		if(sort_every > 0) throw oxDNAException("Polymer interaction is not compatible with particle sorting, aborting");
	}
}

template<typename number, typename number4>
void CUDAPolymerInteraction<number, number4>::cuda_init(number box_side, int N) {
	CUDABaseInteraction<number, number4>::cuda_init(box_side, N);
	PolymerInteraction<number>::init();

	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &N, sizeof(int)) );
	COPY_NUMBER_TO_FLOAT(MD_sqr_rfene, this->_sqr_rfene);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rcut, this->_sqr_rcut);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rep_rcut, this->_Polymer_sqr_rep_rcut);
	COPY_NUMBER_TO_FLOAT(MD_Polymer_lambda, this->_Polymer_lambda);

	if(this->_use_edge) CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_n_forces, &this->_n_forces, sizeof(int)) );
}

template<typename number, typename number4>
void CUDAPolymerInteraction<number, number4>::compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box) {
	CUDASimpleVerletList<number, number4> *_v_lists = dynamic_cast<CUDASimpleVerletList<number, number4> *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) {
			polymer_forces_edge_nonbonded<number, number4>
				<<<(_v_lists->_N_edges - 1)/(this->_launch_cfg.threads_per_block) + 1, this->_launch_cfg.threads_per_block>>>
				(d_poss, this->_d_edge_forces, _v_lists->_d_edge_list, _v_lists->_N_edges, d_box);

			this->_sum_edge_forces(d_forces);

			// potential for removal here
			cudaThreadSynchronize();
			CUT_CHECK_ERROR("forces_second_step error -- after non-bonded");

			polymer_forces_edge_bonded<number, number4>
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss, d_forces, d_bonds);
		}
		else {
			polymer_forces<number, number4>
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss, d_forces, _v_lists->_d_matrix_neighs, _v_lists->_d_number_neighs, d_bonds, d_box);
			CUT_CHECK_ERROR("forces_second_step simple_lists error");
		}
	}
	else {
		CUDANoList<number, number4> *_no_lists = dynamic_cast<CUDANoList<number, number4> *>(lists);

		if(_no_lists != NULL) {
			polymer_forces<number, number4>
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss,  d_forces, d_bonds, d_box);
			CUT_CHECK_ERROR("forces_second_step no_lists error");
		}
	}
}

template class CUDAPolymerInteraction<float, float4>;
template class CUDAPolymerInteraction<double, LR_double4>;
