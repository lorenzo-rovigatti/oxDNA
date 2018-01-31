 /*
 * CUDALevyInteraction.cu
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include "CUDALevyInteraction.h"

#include "CUDA/Lists/CUDASimpleVerletList.h"
#include "CUDA/Lists/CUDANoList.h"
#include "Particles/CustomParticle.h"

/* CUDA constants */
__constant__ bool MD_rigid_model[1];

__constant__ int MD_N[1];
__constant__ int MD_N_centres[1];
__constant__ int MD_patchy_power[1];

__constant__ float MD_sigma[3];
__constant__ float MD_sqr_sigma[3];
__constant__ float MD_sqr_tot_rcut[3];
__constant__ float MD_E_cut[3];

__constant__ float MD_epsilon[1], MD_monomer_epsilon[1];
__constant__ float MD_patch_E_cut[1], MD_patch_monomer_E_cut[1];
__constant__ float MD_sqr_patch_rcut[1];
__constant__ float MD_patch_pow_alpha[1];

__constant__ float MD_fene_K[1];
__constant__ float MD_fene_sqr_r0[1];
__constant__ float MD_lin_k[1];
__constant__ float MD_terminal_lin_k[1];
__constant__ float MD_sqr_rcut[1];

#include "CUDA/cuda_utils/CUDA_lr_common.cuh"

__device__ int get_type_from_btype(int btype) {
	if(btype == LevyInteraction<float>::DIMER_CENTRE) return P_B;
	return P_A;
}

template<typename number, typename number4>
__device__ void _fene(number4 &r, number4 &F) {
	number sqr_r = CUDA_DOT(r, r);
	number energy = -0.5f*MD_fene_K[0]*MD_fene_sqr_r0[0]*logf(1.f - sqr_r/MD_fene_sqr_r0[0]);
	// this number is the module of the force over r, so we don't have to divide the distance
	// vector by its module
	number force_mod = -MD_fene_K[0]*MD_fene_sqr_r0[0]/(MD_fene_sqr_r0[0] - sqr_r);

	F.x -= r.x*force_mod;
	F.y -= r.y*force_mod;
	F.z -= r.z*force_mod;
	F.w += energy;
}

template<typename number, typename number4>
__device__ void _repulsion(number4 &r, int type, number4 &F) {
	number sqr_r = CUDA_DOT(r, r);
	number part = powf(MD_sqr_sigma[type]/sqr_r, MD_patchy_power[0]*0.5f);
	number energy = part - MD_E_cut[type];

	number force_mod = MD_patchy_power[0]*part/sqr_r;

	if(sqr_r > MD_sqr_tot_rcut[type]) energy = force_mod = 0.f;

	F.x -= r.x*force_mod;
	F.y -= r.y*force_mod;
	F.z -= r.z*force_mod;
	F.w += energy;
}

template <typename number, typename number4>
__device__ void _particle_particle_bonded_interaction(number4 &ppos, number4 &qpos, number4 &F, CUDABox<number, number4> *box, bool only_fene=false) {
	int ptype = get_type_from_btype(get_particle_btype<number, number4>(ppos));
	int p_idx = get_particle_index<number, number4>(ppos);
	int qtype = get_type_from_btype(get_particle_btype<number, number4>(qpos));
	int q_idx = get_particle_index<number, number4>(qpos);

	number4 r = box->minimum_image(ppos, qpos);
	if(!only_fene) _repulsion<number, number4>(r, ptype + qtype, F);
	_fene<number, number4>(r, F);
}

template <typename number, typename number4>
__device__ void _particle_particle_interaction(number4 &ppos, number4 &qpos, number4 &ppatch, number4 &qpatch, number4 &F, number4 &T, CUDABox<number, number4> *box) {
	int pbtype = get_particle_btype<number, number4>(ppos);
	int ptype = get_type_from_btype(pbtype);
	int p_idx = get_particle_index<number, number4>(ppos);
	int qbtype = get_particle_btype<number, number4>(qpos);
	int qtype = get_type_from_btype(qbtype);
	int q_idx = get_particle_index<number, number4>(qpos);

	int type = ptype + qtype;
	int btype = pbtype + qbtype;

	number4 r = box->minimum_image(ppos, qpos);

	_repulsion<number, number4>(r, type, F);

	// possible patchy interaction
	if(btype == (LevyInteraction<number>::TETRA_PATCHY + LevyInteraction<number>::DIMER_PATCHY) || btype == (LevyInteraction<number>::TETRA_PATCHY + LevyInteraction<number>::MONOMER) || btype == (LevyInteraction<number>::DIMER_PATCHY + LevyInteraction<number>::MONOMER)) {
		number4 patch_dist = r + qpatch - ppatch;
		number sqr_dist = CUDA_DOT(patch_dist, patch_dist);
		if(sqr_dist < MD_sqr_patch_rcut[0]) {
			number interaction_strength = MD_epsilon[0];
			number E_cut = MD_patch_E_cut[0];
			if(btype != (LevyInteraction<number>::TETRA_PATCHY + LevyInteraction<number>::DIMER_PATCHY)) {
				interaction_strength = MD_monomer_epsilon[0];
				E_cut = MD_patch_monomer_E_cut[0];
			}

			number r8b10 = SQR(SQR(sqr_dist)) / (double) MD_patch_pow_alpha[0];
			number exp_part = -1.001f*interaction_strength*expf(-(number)0.5f*r8b10*sqr_dist);

			F.w += exp_part - E_cut;

			number4 tmp_force = patch_dist*(5.f*exp_part*r8b10);
			F -= tmp_force;
			T -= _cross<number>(ppatch, tmp_force);
		}
	}
}

// forces + second step without lists
template <typename number, typename number4>
__global__ void Levy_forces(number4 *poss, GPU_quat<number> *orientations, number4 *forces, number4 *torques, LR_bonds *bonds, CUDABox<number, number4> *box) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	number4 T = make_number4<number, number4>(0.f, 0.f, 0.f, 0.f);
	number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat<number,number4>(orientations[IND], a1, a2, a3);
	number4 ppatch = a1*0.5f;
	LR_bonds bs = bonds[IND];
	number4 ppos = poss[IND];

	if(bs.n3 != P_INVALID) _particle_particle_bonded_interaction<number, number4>(ppos, poss[bs.n3], F, box);
	if(bs.n5 != P_INVALID) _particle_particle_bonded_interaction<number, number4>(ppos, poss[bs.n5], F, box);

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND && bs.n3 != j && bs.n5 != j) {
			number4 qpos = poss[j];
			get_vectors_from_quat<number,number4>(orientations[j], b1, b2, b3);
			number4 qpatch = b1*0.5f;
			_particle_particle_interaction<number, number4>(ppos, qpos, ppatch, qpatch, F, T, box);
		}
	}

	forces[IND] = F;
	torques[IND] = _vectors_transpose_number4_product(a1, a2, a3, T);
}

// forces + second step with verlet lists
template <typename number, typename number4>
__global__ void Levy_forces(number4 *poss, GPU_quat<number> *orientations, number4 *forces, number4 *torques, int *matrix_neighs, int *number_neighs, LR_bonds *bonds, CUDABox<number, number4> *box) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	number4 T = make_number4<number, number4>(0.f, 0.f, 0.f, 0.f);
	number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat<number,number4>(orientations[IND], a1, a2, a3);
	number4 ppatch = a1*0.5f;
	LR_bonds bs = bonds[IND];
	number4 ppos = poss[IND];

	if(bs.n3 != P_INVALID) _particle_particle_bonded_interaction<number, number4>(ppos, poss[bs.n3], F, box);
	if(bs.n5 != P_INVALID) _particle_particle_bonded_interaction<number, number4>(ppos, poss[bs.n5], F, box);

	const int num_neighs = number_neighs[IND];
	for(int j = 0; j < num_neighs; j++) {
		const int k_index = matrix_neighs[j*MD_N[0] + IND];

		number4 qpos = poss[k_index];
		get_vectors_from_quat<number,number4>(orientations[k_index], b1, b2, b3);
		number4 qpatch = b1*0.5f;
		_particle_particle_interaction<number, number4>(ppos, qpos, ppatch, qpatch, F, T, box);
	}

	forces[IND] = F;
	torques[IND] = _vectors_transpose_number4_product(a1, a2, a3, T);
}

template<typename number, typename number4>
__device__ void _three_body(number4 &ppos, LR_bonds &bs, number4 &F, number4 *poss, number4 *n3_forces, number4 *n5_forces, CUDABox<number, number4> *box) {
	if(bs.n3 == P_INVALID || bs.n5 == P_INVALID) return;

	number4 n3_pos = poss[bs.n3];
	number4 n5_pos = poss[bs.n5];

	int n3_btype = get_particle_btype<number, number4>(n3_pos);
	int n5_btype = get_particle_btype<number, number4>(n5_pos);

	number curr_lin_k = MD_lin_k[0];
	if(n3_btype == LevyInteraction<number>::DIMER_PATCHY || n5_btype == LevyInteraction<number>::DIMER_PATCHY) {
		// if LEVY_rigid != true we don't want the first and last beads to feel any force
		if(!MD_rigid_model[0]) return;
		else curr_lin_k = MD_terminal_lin_k[0];
	}

	number4 dist_pn3 = box->minimum_image(ppos, n3_pos);
	number4 dist_pn5 = box->minimum_image(n5_pos, ppos);

	number sqr_dist_pn3 = CUDA_DOT(dist_pn3, dist_pn3);
	number sqr_dist_pn5 = CUDA_DOT(dist_pn5, dist_pn5);
	number i_pn3_pn5 = 1.f/sqrtf(sqr_dist_pn3*sqr_dist_pn5);
	number cost = CUDA_DOT(dist_pn3, dist_pn5) * i_pn3_pn5;

	number cost_n3 = cost/sqr_dist_pn3;
	number cost_n5 = cost/sqr_dist_pn5;
	number force_mod_n3 = i_pn3_pn5 + cost_n3;
	number force_mod_n5 = i_pn3_pn5 + cost_n5;

	F += dist_pn3*(force_mod_n3*curr_lin_k) - dist_pn5*(force_mod_n5*curr_lin_k);
	// the factor of two comes from the fact that the final energy should be twice the "real" potential energy
	F.w += 2.f*curr_lin_k*(1.f - cost);

	number4 n3_force = dist_pn5*(i_pn3_pn5*curr_lin_k) - dist_pn3*(cost_n3*curr_lin_k);
	number4 n5_force = dist_pn5*(cost_n5*curr_lin_k) - dist_pn3*(i_pn3_pn5*curr_lin_k);
	/*n3_forces[bs.n3] = n3_force;
	  n5_forces[bs.n5] = n5_force;*/
	LR_atomicAddXYZ(n3_forces + bs.n3, n3_force);
	LR_atomicAddXYZ(n5_forces + bs.n5, n5_force);
}

template <typename number, typename number4>
__global__ void three_body_forces(number4 *poss, number4 *forces, number4 *n3_forces, number4 *n5_forces, LR_bonds *bonds, CUDABox<number, number4> *box) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	LR_bonds bs = bonds[IND];
	number4 ppos = poss[IND];

	_three_body<number, number4>(ppos, bs, F, poss, n3_forces, n5_forces, box);

	forces[IND] = F;
}

template <typename number, typename number4>
__global__ void sum_three_body(number4 *forces, number4 *n3_forces, number4 *n5_forces) {
	if(IND >= MD_N[0]) return;

	forces[IND] = forces[IND] + n3_forces[IND] + n5_forces[IND];
}

template <typename number, typename number4>
__global__ void centre_forces(number4 *poss, number4 *forces, number4 *n3_forces, number4 *n5_forces, int *centres, centre_bonds *bonds, CUDABox<number, number4> *box) {
	if(IND >= MD_N_centres[0]) return;

	int idx_centre = centres[IND];
	centre_bonds centre_bonds = bonds[IND];
	number4 pos_centre = poss[idx_centre];
	number4 F = forces[idx_centre];
	LR_bonds bs_centre;

	for(int an = 0; an < CENTRE_N_NEIGHS; an++) {
		int bonded_neigh = centre_bonds.n[an];
		// since bonded neighbours of centre are in the centre's neighbouring list, the LJ interaction between
		// the two, from the point of view of the centre, has been already computed and hence the centre-particle
		// interaction reduces to just the fene
		_particle_particle_bonded_interaction<number, number4>(pos_centre, poss[bonded_neigh], F, box, true);
		bs_centre.n3 = bonded_neigh;

		for(int bn = an + 1; bn < CENTRE_N_NEIGHS; bn++) {
			bs_centre.n5 = centre_bonds.n[bn];
			_three_body<number, number4>(pos_centre, bs_centre, F, poss, n3_forces, n5_forces, box);
		}
	}

	forces[idx_centre] = F;
}

template<typename number, typename number4>
CUDALevyInteraction<number, number4>::CUDALevyInteraction() {
	_d_centres = NULL;
	_d_centre_neighs = NULL;
	_d_n3_forces = _d_n5_forces = NULL;
}

template<typename number, typename number4>
CUDALevyInteraction<number, number4>::~CUDALevyInteraction() {
	if(_d_centres != NULL) {
		CUDA_SAFE_CALL( cudaFree(_d_centres) );
		CUDA_SAFE_CALL( cudaFree(_d_centre_neighs) );
	}

	if(_d_n3_forces != NULL) {
		CUDA_SAFE_CALL( cudaFree(_d_n3_forces) );
		CUDA_SAFE_CALL( cudaFree(_d_n5_forces) );
	}
}

template<typename number, typename number4>
void CUDALevyInteraction<number, number4>::get_settings(input_file &inp) {
	LevyInteraction<number>::get_settings(inp);

	int sort_every;
	if(getInputInt(&inp, "CUDA_sort_every", &sort_every, 0) == KEY_FOUND) {
		if(sort_every > 0) throw oxDNAException("Levy interaction is not compatible with particle sorting, aborting");
	}
}

template<typename number, typename number4>
void CUDALevyInteraction<number, number4>::cuda_init(number box_side, int N) {
	CUDABaseInteraction<number, number4>::cuda_init(box_side, N);
	LevyInteraction<number>::init();

	_setup_centres();

	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number4>(&_d_n3_forces, this->_N*sizeof(number4)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number4>(&_d_n5_forces, this->_N*sizeof(number4)) );
	CUDA_SAFE_CALL( cudaMemset(_d_n3_forces, 0, this->_N*sizeof(number4)) );
	CUDA_SAFE_CALL( cudaMemset(_d_n5_forces, 0, this->_N*sizeof(number4)) );

	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_rigid_model, &this->_rigid_model, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &N, sizeof(int)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_patchy_power, &this->_patchy_power, sizeof(int)) );

	COPY_ARRAY_TO_CONSTANT(MD_sigma, this->_sigma, 3);
	COPY_ARRAY_TO_CONSTANT(MD_sqr_sigma, this->_sqr_sigma, 3);
	COPY_ARRAY_TO_CONSTANT(MD_sqr_tot_rcut, this->_sqr_tot_rcut, 3);
	COPY_ARRAY_TO_CONSTANT(MD_E_cut, this->_E_cut, 3);

	COPY_NUMBER_TO_FLOAT(MD_lin_k, this->_lin_k);
	COPY_NUMBER_TO_FLOAT(MD_terminal_lin_k, this->_terminal_lin_k);
	COPY_NUMBER_TO_FLOAT(MD_fene_K, this->_fene_K);
	COPY_NUMBER_TO_FLOAT(MD_fene_sqr_r0, this->_fene_sqr_r0);
	COPY_NUMBER_TO_FLOAT(MD_sqr_rcut, this->_sqr_rcut);

	COPY_NUMBER_TO_FLOAT(MD_epsilon, this->_epsilon);
	COPY_NUMBER_TO_FLOAT(MD_monomer_epsilon, this->_monomer_epsilon);
	COPY_NUMBER_TO_FLOAT(MD_patch_E_cut, this->_patch_E_cut);
	COPY_NUMBER_TO_FLOAT(MD_patch_monomer_E_cut, this->_patch_monomer_E_cut);
	COPY_NUMBER_TO_FLOAT(MD_sqr_patch_rcut, this->_sqr_patch_rcut);
	COPY_NUMBER_TO_FLOAT(MD_patch_pow_alpha, this->_patch_pow_alpha);
}

template<typename number, typename number4>
void CUDALevyInteraction<number, number4>::_setup_centres() {
	BaseParticle<number> **particles = new BaseParticle<number> *[this->_N];
	LevyInteraction<number>::allocate_particles(particles, this->_N);
	int N_strands;
	LevyInteraction<number>::read_topology(this->_N, &N_strands, particles);

	int N_centres = this->_N_tetramers;
	int *h_centres = new int[N_centres];
	centre_bonds *h_centre_neighs = new centre_bonds[N_centres];

	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<int>(&_d_centres, N_centres*sizeof(int)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<centre_bonds>(&_d_centre_neighs, N_centres*sizeof(centre_bonds)) );

	int rel_idx_centre = 0;
	for(int i = 0; i < this->_N; i++) {
		CustomParticle<number> *p = static_cast<CustomParticle<number> *>(particles[i]);
		if(p->btype == LevyInteraction<number>::TETRA_CENTRE) {
			h_centres[rel_idx_centre] = p->index;

			// now load all the centre_bonds structures by looping over all the bonded neighbours
			int nn = 0;
			for(typename set<CustomParticle<number> *>::iterator it = p->bonded_neighs.begin(); it != p->bonded_neighs.end(); it++) {
				if((*it) != p->n5) {
					h_centre_neighs[rel_idx_centre].n[nn] = (*it)->index;
					nn++;
				}
			}
			for(; nn < CENTRE_N_NEIGHS; nn++) h_centre_neighs[rel_idx_centre].n[nn] = P_INVALID;
			rel_idx_centre++;
		}
	}

	if(rel_idx_centre != N_centres) throw oxDNAException("%d centres found, should have been %d", rel_idx_centre, N_centres);

	CUDA_SAFE_CALL( cudaMemcpy(_d_centres, h_centres, N_centres*sizeof(int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_centre_neighs, h_centre_neighs, N_centres*sizeof(centre_bonds), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N_centres, &N_centres, sizeof(int)) );

	for(int i = 0; i < this->_N; i++) delete particles[i];
	delete[] particles;
	delete[] h_centres;
	delete[] h_centre_neighs;
}

template<typename number, typename number4>
void CUDALevyInteraction<number, number4>::compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box) {
	three_body_forces<number, number4>
		<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
		(d_poss, d_forces, _d_n3_forces, _d_n5_forces, d_bonds, d_box);
	CUT_CHECK_ERROR("three_body_forces error");

	centre_forces<number, number4>
		<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
		(d_poss, d_forces, _d_n3_forces, _d_n5_forces, _d_centres, _d_centre_neighs, d_box);
	CUT_CHECK_ERROR("centre_forces error");

	sum_three_body<number, number4>
		<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
		(d_forces, _d_n3_forces, _d_n5_forces);
	CUT_CHECK_ERROR("sum_three_body error");

	CUDA_SAFE_CALL( cudaMemset(_d_n3_forces, 0, this->_N*sizeof(number4)) );
	CUDA_SAFE_CALL( cudaMemset(_d_n5_forces, 0, this->_N*sizeof(number4)) );

	CUDASimpleVerletList<number, number4> *_v_lists = dynamic_cast<CUDASimpleVerletList<number, number4> *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) throw oxDNAException("use_edge unsupported by LevyInteraction");
		else {
			Levy_forces<number, number4>
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss, d_orientations, d_forces, d_torques, _v_lists->_d_matrix_neighs, _v_lists->_d_number_neighs, d_bonds, d_box);
			CUT_CHECK_ERROR("Levy_forces Verlet Lists error");
		}
	}
	else {
		CUDANoList<number, number4> *_no_lists = dynamic_cast<CUDANoList<number, number4> *>(lists);

		if(_no_lists != NULL) {
			Levy_forces<number, number4>
				<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
				(d_poss, d_orientations, d_forces, d_torques, d_bonds, d_box);
			CUT_CHECK_ERROR("Levy_forces no_lists error");
		}
	}
}

extern "C" IBaseInteraction<float> *make_CUDALevyInteraction_float() {
	return new CUDALevyInteraction<float, float4>();
}

extern "C" IBaseInteraction<double> *make_CUDALevyInteraction_double() {
	return new CUDALevyInteraction<double, LR_double4>();
}

template class CUDALevyInteraction<float, float4>;
template class CUDALevyInteraction<double, LR_double4>;
