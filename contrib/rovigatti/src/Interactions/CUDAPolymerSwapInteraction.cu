/*
 * CUDAPolymerSwapInteraction.cu
 *
 *  Created on: 24/mar/2020
 *      Author: lorenzo
 */

#include "CUDAPolymerSwapInteraction.h"

#include "Particles/CustomParticle.h"
#include "CUDA/Lists/CUDASimpleVerletList.h"
#include "CUDA/Lists/CUDANoList.h"

/* System constants */
__constant__ int MD_N[1];
__constant__ int MD_n[1];
__constant__ float MD_sqr_rep_rcut[3];
__constant__ float MD_sqr_rfene[3];
__constant__ float MD_Kfene[3];
__constant__ float MD_WCA_sigma[3];
__constant__ float MD_sqr_rcut[1];
__constant__ float MD_alpha[1];
__constant__ float MD_beta[1];
__constant__ float MD_gamma[1];

__constant__ float MD_sqr_3b_rcut[1];
__constant__ float MD_3b_sigma[1];
__constant__ float MD_3b_rcut[1];
__constant__ float MD_3b_A_part[1];
__constant__ float MD_3b_B_part[1];

#include "CUDA/cuda_utils/CUDA_lr_common.cuh"

__device__ void _WCA(c_number4 &ppos, c_number4 &qpos, int int_type, c_number4 &F, CUDABox *box) {
	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r);

	c_number energy = 0.f;
	// this c_number is the module of the force over r, so we don't have to divide the distance vector by its module
	c_number force_mod = 0.f;

	if(sqr_r < MD_sqr_rep_rcut[int_type]) {
		c_number part = 1.f;
		c_number ir2_scaled = SQR(MD_WCA_sigma[int_type]) / sqr_r;
		for(int i = 0; i < MD_n[0] / 2; i++) {
			part *= ir2_scaled;
		}
		energy += 4.f * part * (part - 1.f) + 1.f - MD_alpha[0];
		force_mod += 4.f * MD_n[0] * part * (2.f * part - 1.f) / sqr_r;
	}
	/*else {
		energy += 0.5f * MD_alpha[0] * (cosf(MD_gamma[0] * sqr_r + MD_beta[0]) - 1.f);
		force_mod += MD_alpha[0] * MD_gamma[0] * sinf(MD_gamma[0] * sqr_r + MD_beta[0]);
	}*/

	if(sqr_r > MD_sqr_rcut[0]) {
		energy = force_mod = (c_number) 0.f;
	}

	F.x -= r.x * force_mod;
	F.y -= r.y * force_mod;
	F.z -= r.z * force_mod;
	F.w += energy;
}

__device__ void _sticky(c_number4 &ppos, c_number4 &qpos, c_number4 &F, CUDABox *box) {
	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r);

	c_number energy = 0.f;
	// this c_number is the module of the force over r, so we don't have to divide the distance vector by its module
	c_number force_mod = 0.f;

	if(sqr_r < MD_sqr_3b_rcut[0]) {
		c_number r_mod = sqrtf(sqr_r);
		c_number exp_part = expf(MD_3b_sigma[0] / (r_mod - MD_3b_rcut[0]));
		c_number tmp_energy = MD_3b_A_part[0] * exp_part * (MD_3b_B_part[0] / SQR(sqr_r) - 1.);

		energy += tmp_energy;

		force_mod = MD_3b_A_part[0] * exp_part * (4. * MD_3b_B_part[0] / (SQR(sqr_r) * r_mod)) + MD_3b_sigma[0] * tmp_energy / SQR(r_mod - MD_3b_rcut[0]);

		/*
		PSBond p_bond(q, _computed_r, r_mod, tmp_energy);
		PSBond q_bond(p, -_computed_r, r_mod, tmp_energy);

		if(!no_three_body) {
			energy += _three_body(p, p_bond, update_forces);
			energy += _three_body(q, q_bond, update_forces);

			_bonds[p->index].insert(p_bond);
			_bonds[q->index].insert(q_bond);
		}*/
	}

	if(sqr_r > MD_sqr_3b_rcut[0]) {
		energy = force_mod = (c_number) 0.f;
	}

	F.x -= r.x * force_mod;
	F.y -= r.y * force_mod;
	F.z -= r.z * force_mod;
	F.w += energy;
}

__device__ void _FENE(c_number4 &ppos, c_number4 &qpos, int int_type, c_number4 &F, CUDABox *box) {
	c_number sqr_rfene = MD_sqr_rfene[int_type];
	c_number Kfene = MD_Kfene[int_type];

	c_number4 r = box->minimum_image(ppos, qpos);
	c_number sqr_r = CUDA_DOT(r, r);

	c_number energy = -Kfene * sqr_rfene * logf(1.f - sqr_r / sqr_rfene);
	// this c_number is the module of the force over r, so we don't have to divide the distance vector by its module
	c_number force_mod = -2.f * Kfene * sqr_rfene / (sqr_rfene - sqr_r);

	F.x -= r.x * force_mod;
	F.y -= r.y * force_mod;
	F.z -= r.z * force_mod;
	F.w += energy;
}

// this kernel is called "some bonded forces" because it computes the FENE potential between all bonded pair, but the WCA potential
// only between monomer-sticky pairs. The monomer-monomer bonded WCA interaction is handled by the ps_forces function
__global__ void ps_some_bonded_forces(c_number4 *poss, c_number4 *forces, int *bonded_neighs, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 ppos = poss[IND];
	int ptype = get_particle_btype(ppos);

	// the first vale of each column is just the number of bonded neighbours
	int n_bonded_neighs = bonded_neighs[IND];

	for(int i = 1; i <= n_bonded_neighs; i++) {
		int n_idx = bonded_neighs[MD_N[0] * i + IND];
		c_number4 qpos = poss[n_idx];
		int qtype = get_particle_btype(qpos);
		int int_type = ptype + qtype;

		_FENE(ppos, qpos, int_type, F, box);
		if(int_type == 1) {
			_WCA(ppos, qpos, int_type, F, box);
		}
	}

	forces[IND] = F;
}

// forces + second step without lists
__global__ void ps_forces(c_number4 *poss, c_number4 *forces, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 ppos = poss[IND];
	int ptype = get_particle_btype(ppos);

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND) {
			c_number4 qpos = poss[j];
			int qtype = get_particle_btype(qpos);
			int int_type = ptype + qtype;

			if(int_type == 0) {
				_WCA(ppos, qpos, int_type, F, box);
			}
			else if(int_type == 2) {
				_sticky(ppos, qpos, F, box);
			}
		}
	}

	forces[IND] = F;
}

// forces + second step with verlet lists
__global__ void ps_forces(c_number4 *poss, c_number4 *forces, int *matrix_neighs, int *c_number_neighs, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 ppos = poss[IND];

	int num_neighs = c_number_neighs[IND];
	int ptype = get_particle_btype(ppos);

	for(int j = 0; j < num_neighs; j++) {
		int k_index = matrix_neighs[j * MD_N[0] + IND];

		c_number4 qpos = poss[k_index];
		int qtype = get_particle_btype(qpos);
		int int_type = ptype + qtype;

		if(int_type == 0) {
			_WCA(ppos, qpos, int_type, F, box);
		}
		else if(int_type == 2) {
			_sticky(ppos, qpos, F, box);
		}
	}

	forces[IND] = F;
}

CUDAPolymerSwapInteraction::CUDAPolymerSwapInteraction() :
				PolymerSwapInteraction() {
	_d_bonded_neighs = NULL;
}

CUDAPolymerSwapInteraction::~CUDAPolymerSwapInteraction() {
	if(_d_bonded_neighs != NULL) {
		CUDA_SAFE_CALL(cudaFree(_d_bonded_neighs));
	}
}

void CUDAPolymerSwapInteraction::get_settings(input_file &inp) {
	PolymerSwapInteraction::get_settings(inp);
}

void CUDAPolymerSwapInteraction::cuda_init(c_number box_side, int N) {
	CUDABaseInteraction::cuda_init(box_side, N);
	PolymerSwapInteraction::init();

	std::vector<BaseParticle *> particles(_N);
	PolymerSwapInteraction::allocate_particles(particles);
	int tmp_N_strands;
	PolymerSwapInteraction::read_topology(&tmp_N_strands, particles);

	int max_n_neighs = 5;
	int n_elems = (max_n_neighs + 1) * _N;
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&_d_bonded_neighs, n_elems * sizeof(int)));
	std::vector<int> h_bonded_neighs(n_elems);

	for(int i = 0; i < _N; i++) {
		CustomParticle *p = static_cast<CustomParticle *>(particles[i]);
		// start from 1, since the first element will contain the number of bonds
		int nb = 1;
		for(auto q : p->bonded_neighs) {
			if(nb > max_n_neighs) {
				throw oxDNAException("CUDAPolymerSwapInteraction: particle %d has more than %d bonded neighbours", p->index, max_n_neighs);
			}
			h_bonded_neighs[_N * nb + i] = q->index;
			nb++;
		}
		h_bonded_neighs[i] = nb - 1;
	}

	CUDA_SAFE_CALL(cudaMemcpy(_d_bonded_neighs, h_bonded_neighs.data(), n_elems * sizeof(int), cudaMemcpyHostToDevice));
	for(auto particle: particles) {
		delete particle;
	}

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_n, &_PS_n, sizeof(int)));
	COPY_ARRAY_TO_CONSTANT(MD_sqr_rep_rcut, _PS_sqr_rep_rcut, 3)
	COPY_ARRAY_TO_CONSTANT(MD_sqr_rfene, _sqr_rfene, 3);
	COPY_ARRAY_TO_CONSTANT(MD_Kfene, _Kfene, 3);
	COPY_ARRAY_TO_CONSTANT(MD_WCA_sigma, _WCA_sigma, 3);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rcut, _sqr_rcut);
	COPY_NUMBER_TO_FLOAT(MD_alpha, _PS_alpha);
	COPY_NUMBER_TO_FLOAT(MD_beta, _PS_beta);
	COPY_NUMBER_TO_FLOAT(MD_gamma, _PS_gamma);
	COPY_NUMBER_TO_FLOAT(MD_sqr_3b_rcut, _sqr_3b_rcut);
	COPY_NUMBER_TO_FLOAT(MD_3b_sigma, _3b_sigma);
	COPY_NUMBER_TO_FLOAT(MD_3b_rcut, _3b_rcut);
	COPY_NUMBER_TO_FLOAT(MD_3b_A_part, _3b_A_part);
	COPY_NUMBER_TO_FLOAT(MD_3b_B_part, _3b_B_part);
}

void CUDAPolymerSwapInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
	ps_some_bonded_forces
		<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
		(d_poss, d_forces, _d_bonded_neighs, d_box);
	CUT_CHECK_ERROR("forces_second_step MG simple_lists error");

	CUDASimpleVerletList *_v_lists = dynamic_cast<CUDASimpleVerletList *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) throw oxDNAException("use_edge unsupported by CPMixtureInteraction");

		ps_forces
			<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
			(d_poss, d_forces, _v_lists->_d_matrix_neighs, _v_lists->_d_c_number_neighs, d_box);
		CUT_CHECK_ERROR("forces_second_step MG simple_lists error");
	}

	CUDANoList *_no_lists = dynamic_cast<CUDANoList *>(lists);
	if(_no_lists != NULL) {
		ps_forces
			<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
			(d_poss, d_forces, d_box);
		CUT_CHECK_ERROR("forces_second_step MG no_lists error");
	}

//	number energy = GpuUtils::sum_c_number4_to_double_on_GPU(d_forces, _N);
//	printf("U: %lf\n", energy / _N / 2.);
}
