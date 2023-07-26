/*
 * CUDATEPInteraction.cu
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#include "CUDATEPInteraction.h"

/* System constants */
__constant__ int MD_N[1];
__constant__ float MD_sqr_rcut[1];
__constant__ float MD_ka[1];
__constant__ float MD_kb[1];
__constant__ float MD_kt[1];

// parameters for the boundary twisting
__constant__ float MD_twist_boundary_stiff[1];
__constant__ float MD_o_modulus[2];

// parameters of the modified twisting terms
__constant__ float MD_twist_a[1];
__constant__ float MD_twist_b[1];
// FENE parameters
__constant__ float MD_TEP_FENE_DELTA[1];
__constant__ float MD_TEP_FENE_DELTA2[1];
__constant__ float MD_TEP_FENE_EPS[1];
__constant__ float MD_TEP_FENE_R0[1];
// LJ parameters
__constant__ float MD_TEP_EXCL_EPS_BONDED[1];
__constant__ float MD_TEP_EXCL_EPS_NONBONDED[1];
__constant__ float MD_TEP_EXCL_S2[1];
__constant__ float MD_TEP_EXCL_R2[1];
__constant__ float MD_TEP_EXCL_B2[1];
__constant__ float MD_TEP_EXCL_RC2[1];
__constant__ float MD_TEP_spring_offset[1];

#include "../cuda_utils/CUDA_lr_common.cuh"

__device__ void _repulsive_lj2(c_number prefactor, c_number4 &r, c_number4 &F, c_number sigma, c_number rstar, c_number b, c_number rc) {
	c_number rnorm = CUDA_DOT(r, r);
	if(rnorm < SQR(rc)) {
		if(rnorm > SQR(rstar)) {
			c_number rmod = sqrtf(rnorm);
			c_number rrc = rmod - rc;
			c_number part = prefactor * b * rrc;
			F.x += r.x * (2.f * part / rmod);
			F.y += r.y * (2.f * part / rmod);
			F.z += r.z * (2.f * part / rmod);
			F.w += part * rrc;
		}
		else {
			c_number tmp = SQR(sigma) / rnorm;
			c_number lj_part = tmp * tmp * tmp;
			// the additive term was added by me to mimick Davide's implementation
			F.x += r.x * (24.f * prefactor * (lj_part - 2.f * SQR(lj_part)) / rnorm);
			F.y += r.y * (24.f * prefactor * (lj_part - 2.f * SQR(lj_part)) / rnorm);
			F.z += r.z * (24.f * prefactor * (lj_part - 2.f * SQR(lj_part)) / rnorm);
			F.w += 4.f * prefactor * (SQR(lj_part) - lj_part) + prefactor;
		}
	}
}

__device__ void _bonded_double_bending(c_number4 &up, c_number4 &uq, c_number kb1, c_number kb2, c_number tu, c_number tk, c_number4 &F, c_number4 &T) {
	c_number4 torque;

	// end of added block
	c_number cosine = CUDA_DOT(up, uq);

	c_number gu = cosf(tu);
	c_number angle = LRACOS(cosine);

	c_number A = (kb2 * sinf(tk) - kb1 * sinf(tu)) / (tk - tu);
	c_number B = kb1 * sinf(tu);
	c_number C = kb1 * (1.f - gu) - (A * SQR(tu) / 2.f + tu * (B - tu * A));
	c_number D = cosf(tk) + (A * SQR(tk) * 0.5f + tk * (B - tu * A) + C) / kb2;

	c_number g1 = (angle - tu) * A + B;
	c_number g_ = A * SQR(angle) / 2.f + angle * (B - tu * A);
	c_number g = g_ + C;

	// unkinked bending regime
	if(angle < tu) {
		torque = MD_kb[0] * kb1 * _cross(up, uq);
		F.w += MD_kb[0] * (1.f - cosine) * kb1;
	}
	// intermediate regime
	else if(angle < tk) {
		c_number4 sin_vector = _cross(up, uq);
		c_number sin_module = _module(sin_vector);
		sin_vector.x /= sin_module;
		sin_vector.y /= sin_module;
		sin_vector.z /= sin_module;

		torque = MD_kb[0] * g1 * sin_vector;
		F.w += MD_kb[0] * g;
	}
	// kinked bending regime - same as unkinked, but with an additive term to the energy and with a different bending.
	else {
		torque = MD_kb[0] * kb2 * _cross(up, uq);
		F.w += MD_kb[0] * (D - cosine) * kb2;
	}

	T.x += torque.x;
	T.y += torque.y;
	T.z += torque.z;
}

template<typename c_number, typename c_number4, bool torque_on_p>
__device__ void _bonded_particle_particle(c_number4 &n3_pos, c_number4 &up, c_number4 &fp, c_number4 &vp, c_number4 &n5_pos, c_number4 &uq, c_number4 &fq, c_number4 &vq, c_number4 &F, c_number4 &T, CUDABox *box, c_number kb1, c_number kb2, c_number tu, c_number tk, c_number kt_pref, bool alignment_only = false) {
	c_number4 r = box->minimum_image(n3_pos, n5_pos);
	c_number rmod = sqrtf(CUDA_DOT(r, r));

	if(!alignment_only) {
		// excluded volume
		_repulsive_lj2(MD_TEP_EXCL_EPS_BONDED[0], r, F, MD_TEP_EXCL_S2[0], MD_TEP_EXCL_R2[0], MD_TEP_EXCL_B2[0], MD_TEP_EXCL_RC2[0]);
		F.w += MD_TEP_spring_offset[0];

		// spring
		c_number r0 = rmod - MD_TEP_FENE_R0[0];
		c_number force_module = (MD_TEP_FENE_EPS[0] * r0 / (MD_TEP_FENE_DELTA2[0] - SQR(r0))) / rmod;
		F.x += r.x * force_module;
		F.y += r.y * force_module;
		F.z += r.z * force_module;
		F.w += -MD_TEP_FENE_EPS[0] * 0.5f * logf(1.f - SQR(r0) / MD_TEP_FENE_DELTA2[0]);

		// bonded twist
		c_number M = CUDA_DOT(fp, fq) + CUDA_DOT(vp, vq);
		c_number L = 1.f + CUDA_DOT(up, uq);
		c_number cos_alpha_plus_gamma = M / L;
		c_number energy = 0.;

		if(-cos_alpha_plus_gamma >= MD_twist_b[0]) printf("We have an issue with the twist angle (%lf) in thread %d\n", cos_alpha_plus_gamma, IND);
		else {
			c_number4 torque = (MD_kt[0] / L) * (_cross(fp, fq) + _cross(vp, vq) - cos_alpha_plus_gamma * _cross(up, uq)) * kt_pref;
			if(-cos_alpha_plus_gamma <= MD_twist_a[0]) {
				energy += MD_kt[0] * (1.f - cos_alpha_plus_gamma);
			}
			else {
				c_number A = SQR(MD_twist_b[0] - MD_twist_a[0]) * CUB(MD_twist_b[0] - MD_twist_a[0]) * 0.25f;	//A = (b - a)^5/4
				c_number C = 0.25f * (MD_twist_a[0] - MD_twist_b[0]) + MD_twist_a[0];	// C = (a - b)/4 + a
				torque *= 4.f * A / (SQR(MD_twist_b[0] + cos_alpha_plus_gamma) * CUB(MD_twist_b[0] + cos_alpha_plus_gamma));
				energy += MD_kt[0] * (1.f + A / (CUB(-cos_alpha_plus_gamma - MD_twist_b[0])) + C);
			}

			T.x += torque.x;
			T.y += torque.y;
			T.z += torque.z;
		}

		F.w += kt_pref * energy;

		// bonded bending
		_bonded_double_bending(up, uq, kb1, kb2, tu, tk, F, T);
		// this is the old "single" bending
		/*c_number4 torque = MD_kb[0]*kb1*_cross(up, uq);
		 T.x += torque.x;
		 T.y += torque.y;
		 T.z += torque.z;
		 F.w += MD_kb[0]*(1.f - CUDA_DOT(up, uq))*kb1;*/
	}

	// bonded alignment
	// here we have to take into account the intrinsic directionality of the model (i.e. forces and torques have different expressions if we go from 3' to 5' or vice versa)
	bool keep_r = ((torque_on_p && !alignment_only) || (!torque_on_p && alignment_only));
	c_number4 ba_tp = (keep_r) ? r : -r;
	c_number4 ba_up = (torque_on_p) ? up : uq;
	c_number4 force = MD_ka[0] * (ba_up - ba_tp * CUDA_DOT(ba_up, ba_tp) / SQR(rmod)) / rmod;
	if(keep_r) force = -force;
	F.x += force.x;
	F.y += force.y;
	F.z += force.z;
	F.w += MD_ka[0] * (1.f - CUDA_DOT(ba_up, ba_tp) / rmod);
	if(torque_on_p) {
		c_number4 torque = -(MD_ka[0] * _cross(ba_tp, ba_up)) / rmod;
		T.x += torque.x;
		T.y += torque.y;
		T.z += torque.z;
	}
}

__device__ c_number4 rotate_vector_around_versor(c_number4 &vector, c_number4 &versor, c_number angle) {
	c_number costh = cosf(angle);
	c_number sinth = sinf(angle);
	c_number scalar = CUDA_DOT(vector, versor);
	c_number4 cross = _cross(versor, vector);
	return make_c_number4(versor.x * scalar * (1.f - costh) + vector.x * costh + cross.x * sinth, versor.y * scalar * (1.f - costh) + vector.y * costh + cross.y * sinth, versor.z * scalar * (1.f - costh) + vector.z * costh + cross.z * sinth, 0.f);
}

__device__ void _twist_boundary_particle(c_number4 &v1, c_number4& v3, c_number4 &T, c_number4 &o_vect, c_number4 &w_vect, llint step, c_number o_modulus) {
	if(MD_twist_boundary_stiff[0] == 0.f) return;
	c_number4 wt = rotate_vector_around_versor(w_vect, o_vect, step * o_modulus);

	c_number4 torque_w = MD_twist_boundary_stiff[0] * (_cross(v3, wt));
	c_number4 torque_o = MD_twist_boundary_stiff[0] * (_cross(v1, o_vect));
	T += torque_w + torque_o;
}

__global__ void TEP_forces(c_number4 *poss, GPU_quat *orientations, c_number4 *forces, c_number4 *torques, int *matrix_neighs, int *number_neighs, LR_bonds *bonds, CUDABox *box, c_number *kb1_pref, c_number *kb2_pref, c_number *xk_bending, c_number *xu_bending, c_number *kt_pref, c_number4 *o_vects, c_number4 *w_vects, llint step) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = make_c_number4(0, 0, 0, 0);
	c_number4 ppos = poss[IND];
	LR_bonds pbonds = bonds[IND];

	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[IND], a1, a2, a3);

	if(pbonds.n3 != P_INVALID) {
		int n3_index = pbonds.n3;
		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[pbonds.n3], b1, b2, b3);

		c_number kb1 = kb1_pref[n3_index];
		c_number kb2 = kb2_pref[n3_index];
		c_number xk = xk_bending[n3_index];
		c_number xu = xu_bending[n3_index];
		c_number kt = kt_pref[n3_index];
		_bonded_particle_particle<c_number, c_number4, false>(ppos, a1, a2, a3, poss[n3_index], b1, b2, b3, F, T, box, kb1, kb2, xk, xu, kt);
		// if this bead is the last of the strand
		if(pbonds.n5 == P_INVALID) _bonded_particle_particle<c_number, c_number4, true>(ppos, a1, a2, a3, poss[pbonds.n3], b1, b2, b3, F, T, box, kb1, kb2, xk, xu, kt, true);
	}
	else _twist_boundary_particle(a1, a3, T, o_vects[0], w_vects[0], step, MD_o_modulus[0]);

	if(pbonds.n5 != P_INVALID) {
		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[pbonds.n5], b1, b2, b3);

		c_number kb1 = kb1_pref[IND];
		c_number kb2 = kb2_pref[IND];
		c_number xk = xk_bending[IND];
		c_number xu = xu_bending[IND];
		c_number kt = kt_pref[IND];
		_bonded_particle_particle<c_number, c_number4, true>(ppos, a1, a2, a3, poss[pbonds.n5], b1, b2, b3, F, T, box, kb1, kb2, xk, xu, kt);
		// we have to check whether n5 is the last of the strand
		LR_bonds qbonds = bonds[pbonds.n5];
		if(qbonds.n5 == P_INVALID) {
			_bonded_particle_particle<c_number, c_number4, false>(ppos, a1, a2, a3, poss[pbonds.n5], b1, b2, b3, F, T, box, kb1, kb2, xk, xu, kt, true);
			_twist_boundary_particle(a1, a3, T, o_vects[1], w_vects[1], step, MD_o_modulus[1]);
		}
	}

	int num_neighs = NUMBER_NEIGHBOURS(IND, number_neighs);
	for(int j = 0; j < num_neighs; j++) {
		int k_index = NEXT_NEIGHBOUR(IND, j, matrix_neighs);

		if(k_index != IND) {
			c_number4 qpos = poss[k_index];
			c_number4 r = box->minimum_image(ppos, qpos);
			_repulsive_lj2(MD_TEP_EXCL_EPS_NONBONDED[0], r, F, MD_TEP_EXCL_S2[0], MD_TEP_EXCL_R2[0], MD_TEP_EXCL_B2[0], MD_TEP_EXCL_RC2[0]);
		}
	}

	F.w *= (c_number) 0.5f;
	forces[IND] = F;

	T = _vectors_transpose_c_number4_product(a1, a2, a3, T);
	torques[IND] = T;
}

#include "../Lists/CUDASimpleVerletList.h"
#include "../Lists/CUDANoList.h"

CUDATEPInteraction::CUDATEPInteraction() {
	_d_kb1_pref = _d_kb2_pref = _d_kt_pref = _d_xk_bending = _d_xu_bending = NULL;
	_d_o_vects = _d_w_vects = NULL;
	_steps = 0;
}

CUDATEPInteraction::~CUDATEPInteraction() {
	if(_d_kb1_pref != NULL) CUDA_SAFE_CALL(cudaFree(_d_kb1_pref));
	if(_d_kb2_pref != NULL) CUDA_SAFE_CALL(cudaFree(_d_kb2_pref));
	if(_d_xk_bending != NULL) CUDA_SAFE_CALL(cudaFree(_d_xk_bending));
	if(_d_xu_bending != NULL) CUDA_SAFE_CALL(cudaFree(_d_xu_bending));
	if(_d_kt_pref != NULL) CUDA_SAFE_CALL(cudaFree(_d_kt_pref));
	if(_d_o_vects != NULL) CUDA_SAFE_CALL(cudaFree(_d_o_vects));
	if(_d_w_vects != NULL) CUDA_SAFE_CALL(cudaFree(_d_w_vects));
}

void CUDATEPInteraction::get_settings(input_file &inp) {
	TEPInteraction::get_settings(inp);

	int sort_every;
	if(getInputInt(&inp, "CUDA_sort_every", &sort_every, 0) == KEY_FOUND) {
		if(sort_every > 0) throw oxDNAException("TEP interaction is not compatible with particle sorting, aborting");
	}

	if(this->_prefer_harmonic_over_fene) throw oxDNAException("The 'prefer_harmonic_over_fene' option is not compatible with CUDA");
}

void CUDATEPInteraction::cuda_init(int N) {
	CUDABaseInteraction::cuda_init(N);
	TEPInteraction::init();

	std::vector<BaseParticle *> particles(N);
	int my_N_strands;
	TEPInteraction::read_topology(&my_N_strands, particles);

	size_t k_size = N * sizeof(c_number);
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_kb1_pref, k_size));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_kb2_pref, k_size));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_xk_bending, k_size));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_xu_bending, k_size));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc(&_d_kt_pref, k_size));

	CUDA_SAFE_CALL(cudaMemcpy(_d_kb1_pref, this->_kb1_pref, k_size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(_d_kb2_pref, this->_kb2_pref, k_size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(_d_xk_bending, this->_xk_bending, k_size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(_d_xu_bending, this->_xu_bending, k_size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(_d_kt_pref, this->_kt_pref, k_size, cudaMemcpyHostToDevice));

	for(int i = 0; i < N; i++) {
		delete particles[i];
	}

	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<c_number4>(&_d_o_vects, 2 * sizeof(c_number4)));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<c_number4>(&_d_w_vects, 2 * sizeof(c_number4)));
	c_number4 h_o_vects[2] = { { (c_number) _o1.x, (c_number) _o1.y, (c_number) _o1.z, (float) 0. }, { (c_number) _o2.x, (c_number) _o2.y, (c_number) _o2.z, (float) 0. } };
	c_number4 h_w_vects[2] = { { (c_number) _w1.x, (c_number) _w1.y, (c_number) _w1.z, (float) 0. }, { (c_number) _w2.x, (c_number) _w2.y, (c_number) _w2.z, (float) 0. } };

	CUDA_SAFE_CALL(cudaMemcpy(_d_o_vects, h_o_vects, 2 * sizeof(c_number4), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(_d_w_vects, h_w_vects, 2 * sizeof(c_number4), cudaMemcpyHostToDevice));
	COPY_NUMBER_TO_FLOAT(MD_twist_boundary_stiff, this->_twist_boundary_stiff);
	float o_modulus[2] = { (float) this->_o1_modulus, (float) this->_o2_modulus };
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_o_modulus, o_modulus, 2 * sizeof(float)));

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));

	COPY_NUMBER_TO_FLOAT(MD_sqr_rcut, this->_sqr_rcut);
	COPY_NUMBER_TO_FLOAT(MD_ka, this->_ka);
	COPY_NUMBER_TO_FLOAT(MD_kb, this->_kb);
	COPY_NUMBER_TO_FLOAT(MD_kt, this->_kt);
	COPY_NUMBER_TO_FLOAT(MD_twist_a, this->_twist_a);
	COPY_NUMBER_TO_FLOAT(MD_twist_b, this->_twist_b);
	COPY_NUMBER_TO_FLOAT(MD_TEP_FENE_DELTA, this->_TEP_FENE_DELTA);
	COPY_NUMBER_TO_FLOAT(MD_TEP_FENE_DELTA2, this->_TEP_FENE_DELTA2);
	COPY_NUMBER_TO_FLOAT(MD_TEP_FENE_EPS, this->_TEP_FENE_EPS);
	COPY_NUMBER_TO_FLOAT(MD_TEP_FENE_R0, this->_TEP_FENE_R0);
	COPY_NUMBER_TO_FLOAT(MD_TEP_EXCL_EPS_BONDED, this->_TEP_EXCL_EPS_BONDED);
	COPY_NUMBER_TO_FLOAT(MD_TEP_EXCL_EPS_NONBONDED, this->_TEP_EXCL_EPS_NONBONDED);
	COPY_NUMBER_TO_FLOAT(MD_TEP_EXCL_S2, this->_TEP_EXCL_S2);
	COPY_NUMBER_TO_FLOAT(MD_TEP_EXCL_R2, this->_TEP_EXCL_R2);
	COPY_NUMBER_TO_FLOAT(MD_TEP_EXCL_B2, this->_TEP_EXCL_B2);
	COPY_NUMBER_TO_FLOAT(MD_TEP_EXCL_RC2, this->_TEP_EXCL_RC2);
	COPY_NUMBER_TO_FLOAT(MD_TEP_spring_offset, this->_TEP_spring_offset);
}

void CUDATEPInteraction::compute_forces(CUDABaseList*lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox*d_box) {
	_steps++;

	update_increment(_steps);
	_time_var = this->update_time_variable(_time_var);

	TEP_forces
		<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
		(d_poss, d_orientations, d_forces, d_torques, lists->d_matrix_neighs, lists->d_number_neighs, d_bonds, d_box, _d_kb1_pref, _d_kb2_pref, _d_xk_bending, _d_xu_bending, _d_kt_pref, _d_o_vects, _d_w_vects, this->_time_var);
	CUT_CHECK_ERROR("forces_second_step TEP simple_lists error");
}
