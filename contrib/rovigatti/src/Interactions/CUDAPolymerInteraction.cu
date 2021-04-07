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

__device__ void _nonbonded(c_number4 &r, int int_type, c_number4 &F) {
	c_number sqr_r = CUDA_DOT(r, r);

	c_number energy = 0.f;
	// this c_number is the module of the force over r, so we don't have to divide the distance
	// vector for its module
	c_number force_mod = 0.f;

	if(sqr_r < MD_sqr_rep_rcut[0]) {
		c_number part = powf(1.f / sqr_r, 3.f);
		energy += 4.f * part * (part - 1.f) + 1.f;
		force_mod += 24.f * part * (2.f * part - 1.f) / sqr_r;
	}

	if(int_type == 2) {
		if(sqr_r < MD_sqr_rep_rcut[0]) energy -= MD_Polymer_lambda[0];
		else {
			c_number part = powf(1.f / sqr_r, 3.f);
			energy += 24.f * part * (part - 1.f);
			force_mod += 24.f * MD_Polymer_lambda[0] * part * (2 * part - 1.f) / sqr_r;
		}
	}

	if(sqr_r > MD_sqr_rcut[0]) energy = force_mod = (c_number) 0.f;

	F.x -= r.x * force_mod;
	F.y -= r.y * force_mod;
	F.z -= r.z * force_mod;
	F.w += energy;
}

__device__ void _fene(c_number4 &r, c_number4 &F) {
	c_number sqr_r = CUDA_DOT(r, r);

	c_number energy = -15.f * MD_sqr_rfene[0] * logf(1.f - sqr_r / MD_sqr_rfene[0]);

	// this c_number is the module of the force over r, so we don't have to divide the distance
	// vector by its module
	c_number force_mod = -30.f * MD_sqr_rfene[0] / (MD_sqr_rfene[0] - sqr_r);
	F.x -= r.x * force_mod;
	F.y -= r.y * force_mod;
	F.z -= r.z * force_mod;
	F.w += energy;
}

__device__ void _particle_particle_bonded_interaction(c_number4 &ppos, c_number4 &qpos, c_number4 &F) {
	c_number4 r = qpos - ppos;
	// bonded interactions are purely repulsive, so we set the int_type to 0 for all pairs
	_nonbonded(r, 0, F);
	_fene(r, F);
}

__device__ void _particle_particle_interaction(c_number4 &ppos, c_number4 &qpos, c_number4 &F, CUDABox *box) {
	int ptype = get_particle_type(ppos);
	int qtype = get_particle_type(qpos);
	int int_type = ptype + qtype;

	c_number4 r = box->minimum_image(ppos, qpos);
	_nonbonded(r, int_type, F);
}

// forces + second step without lists

__global__ void polymer_forces(c_number4 *poss, c_number4 *forces, LR_bonds *bonds, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	LR_bonds bs = bonds[IND];
	c_number4 ppos = poss[IND];

	if(bs.n3 != P_INVALID) {
		c_number4 qpos = poss[bs.n3];
		_particle_particle_bonded_interaction(ppos, qpos, F);
	}

	if(bs.n5 != P_INVALID) {
		c_number4 qpos = poss[bs.n5];
		_particle_particle_bonded_interaction(ppos, qpos, F);
	}

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND && bs.n3 != j && bs.n5 != j) {
			c_number4 qpos = poss[j];
			_particle_particle_interaction(ppos, qpos, F, box);
		}
	}

	forces[IND] = F;
}

__global__ void polymer_forces_edge_nonbonded(c_number4 *poss, c_number4 *forces, edge_bond *edge_list, int n_edges, CUDABox *box) {
	if(IND >= n_edges) return;

	c_number4 dF = make_c_number4(0, 0, 0, 0);

	edge_bond b = edge_list[IND];

	// get info for particle 1
	c_number4 ppos = poss[b.from];

	// get info for particle 2
	c_number4 qpos = poss[b.to];

	_particle_particle_interaction(ppos, qpos, dF, box);

	dF.w *= (c_number) 0.5f;

	int from_index = MD_N[0] * (IND % MD_n_forces[0]) + b.from;
	//int from_index = MD_N[0]*(b.n_from % MD_n_forces[0]) + b.from;
	if((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (c_number) 0.f) LR_atomicAddXYZ(&(forces[from_index]), dF);

	// Allen Eq. 6 pag 3:
	c_number4 dr = box->minimum_image(ppos, qpos); // returns qpos-ppos
	c_number4 crx = _cross(dr, dF);

	dF.x = -dF.x;
	dF.y = -dF.y;
	dF.z = -dF.z;

	int to_index = MD_N[0] * (IND % MD_n_forces[0]) + b.to;
	if((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (c_number) 0.f) LR_atomicAddXYZ(&(forces[to_index]), dF);
}

// bonded interactions for edge-based approach

__global__ void polymer_forces_edge_bonded(c_number4 *poss, c_number4 *forces, LR_bonds *bonds) {
	if(IND >= MD_N[0]) return;

	c_number4 F0;

	F0.x = forces[IND].x;
	F0.y = forces[IND].y;
	F0.z = forces[IND].z;
	F0.w = forces[IND].w;

	c_number4 dF = make_c_number4(0, 0, 0, 0);
	c_number4 ppos = poss[IND];
	LR_bonds bs = bonds[IND];

	if(bs.n3 != P_INVALID) {
		c_number4 qpos = poss[bs.n3];
		_particle_particle_bonded_interaction(ppos, qpos, dF);
	}
	if(bs.n5 != P_INVALID) {
		c_number4 qpos = poss[bs.n5];
		_particle_particle_bonded_interaction(ppos, qpos, dF);
	}

	forces[IND] = (dF + F0);
}

// forces + second step with verlet lists

__global__ void polymer_forces(c_number4 *poss, c_number4 *forces, int *matrix_neighs, int *c_number_neighs, LR_bonds *bonds, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 ppos = poss[IND];
	LR_bonds bs = bonds[IND];

	if(bs.n3 != P_INVALID) {
		c_number4 qpos = poss[bs.n3];
		_particle_particle_bonded_interaction(ppos, qpos, F);
	}
	if(bs.n5 != P_INVALID) {
		c_number4 qpos = poss[bs.n5];
		_particle_particle_bonded_interaction(ppos, qpos, F);
	}

	const int num_neighs = c_number_neighs[IND];
	for(int j = 0; j < num_neighs; j++) {
		const int k_index = matrix_neighs[j * MD_N[0] + IND];

		c_number4 qpos = poss[k_index];
		_particle_particle_interaction(ppos, qpos, F, box);
	}

	forces[IND] = F;
}

/* END CUDA */

CUDAPolymerInteraction::CUDAPolymerInteraction() {

}

CUDAPolymerInteraction::~CUDAPolymerInteraction() {

}

void CUDAPolymerInteraction::get_settings(input_file &inp) {
	PolymerInteraction::get_settings(inp);

	int sort_every;
	if(getInputInt(&inp, "CUDA_sort_every", &sort_every, 0) == KEY_FOUND) {
		if(sort_every > 0) throw oxDNAException("Polymer interaction is not compatible with particle sorting, aborting");
	}
}

void CUDAPolymerInteraction::cuda_init(c_number box_side, int N) {
	CUDABaseInteraction::cuda_init(box_side, N);
	PolymerInteraction::init();

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));
	COPY_NUMBER_TO_FLOAT(MD_sqr_rfene, this->_sqr_rfene);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rcut, this->_sqr_rcut);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rep_rcut, this->_Polymer_sqr_rep_rcut);
	COPY_NUMBER_TO_FLOAT(MD_Polymer_lambda, this->_Polymer_lambda);

	if(this->_use_edge) CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_n_forces, &this->_n_forces, sizeof(int)));
}

void CUDAPolymerInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
	CUDASimpleVerletList *_v_lists = dynamic_cast<CUDASimpleVerletList *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) {
			polymer_forces_edge_nonbonded
				<<<(_v_lists->N_edges - 1)/(this->_launch_cfg.threads_per_block) + 1, this->_launch_cfg.threads_per_block>>>
				(d_poss, this->_d_edge_forces, _v_lists->d_edge_list, _v_lists->N_edges, d_box);

			this->_sum_edge_forces(d_forces);

			// potential for removal here
			cudaThreadSynchronize();
			CUT_CHECK_ERROR("forces_second_step error -- after non-bonded");

			polymer_forces_edge_bonded
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss, d_forces, d_bonds);
		}
		else {
			polymer_forces
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss, d_forces, _v_lists->d_matrix_neighs, _v_lists->d_number_neighs, d_bonds, d_box);
			CUT_CHECK_ERROR("forces_second_step simple_lists error");
		}
	}
	else {
		CUDANoList *_no_lists = dynamic_cast<CUDANoList *>(lists);

		if(_no_lists != NULL) {
			polymer_forces
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss,  d_forces, d_bonds, d_box);
			CUT_CHECK_ERROR("forces_second_step no_lists error");
		}
	}
}
