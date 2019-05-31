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

template<typename number, typename number4>
__device__ void _repulsive_lj2(number prefactor, number4 &r, number4 &F, number sigma, number rstar, number b, number rc) {
	number rnorm = CUDA_DOT(r, r);
	if(rnorm < SQR(rc)) {
		if(rnorm > SQR(rstar)) {
			number rmod = sqrtf(rnorm);
			number rrc = rmod - rc;
			number part = prefactor*b*rrc;
			F.x += r.x*(2.f*part/rmod);
			F.y += r.y*(2.f*part/rmod);
			F.z += r.z*(2.f*part/rmod);
			F.w += part*rrc;
		}
		else {
			number tmp = SQR(sigma)/rnorm;
			number lj_part = tmp*tmp*tmp;
			// the additive term was added by me to mimick Davide's implementation
			F.x += r.x*(24.f*prefactor*(lj_part - 2.f*SQR(lj_part))/rnorm);
			F.y += r.y*(24.f*prefactor*(lj_part - 2.f*SQR(lj_part))/rnorm);
			F.z += r.z*(24.f*prefactor*(lj_part - 2.f*SQR(lj_part))/rnorm);
			F.w += 4.f*prefactor*(SQR(lj_part) - lj_part) + prefactor;
		}
	}
}

template<typename number, typename number4>
__device__ void _bonded_double_bending(number4 &up, number4 &uq, number kb1, number kb2, number tu, number tk, number4 &F, number4 &T) {
	number4 torque;

	// end of added block
	number cosine = CUDA_DOT(up, uq);

	number gu = cosf(tu);
	number angle = LRACOS(cosine);

	number A = (kb2*sinf(tk) - kb1*sinf(tu))/(tk - tu);
	number B = kb1*sinf(tu);
	number C = kb1*(1.f - gu) - (A*SQR(tu)/2.f + tu*(B - tu*A));
	number D = cosf(tk) + (A*SQR(tk)*0.5f + tk * (B - tu*A) + C)/kb2;

	number g1 = (angle - tu)*A + B;
	number g_ = A*SQR(angle)/2.f + angle*(B - tu*A);
	number g = g_ + C;

	// unkinked bending regime
	if(angle < tu) {
		torque = MD_kb[0]*kb1*_cross<number, number4>(up, uq);
		F.w += MD_kb[0]*(1.f - cosine)*kb1;
	}
	// intermediate regime
	else if(angle < tk) {
		number4 sin_vector = _cross<number, number4>(up, uq);
		number sin_module = _module<number, number4>(sin_vector);
		sin_vector.x /= sin_module;
		sin_vector.y /= sin_module;
		sin_vector.z /= sin_module;

		torque = MD_kb[0]*g1*sin_vector;
		F.w += MD_kb[0]*g;
	}
	// kinked bending regime - same as unkinked, but with an additive term to the energy and with a different bending.
	else {
		torque = MD_kb[0]*kb2*_cross<number, number4>(up, uq);
		F.w += MD_kb[0]*(D - cosine)*kb2;
	}

	T.x += torque.x;
	T.y += torque.y;
	T.z += torque.z;
}

template<typename number, typename number4, bool torque_on_p>
__device__ void _bonded_particle_particle(number4 &n3_pos, number4 &up, number4 &fp, number4 &vp,
										number4 &n5_pos, number4 &uq, number4 &fq, number4 &vq,
										number4 &F, number4 &T, CUDABox<number, number4> *box,
								        number kb1, number kb2, number tu, number tk, number kt_pref,
										bool alignment_only=false) {
	number4 r = box->minimum_image(n3_pos, n5_pos);
	number rmod = sqrtf(CUDA_DOT(r, r));

	if(!alignment_only) {
		// excluded volume
		_repulsive_lj2(MD_TEP_EXCL_EPS_BONDED[0], r, F, MD_TEP_EXCL_S2[0], MD_TEP_EXCL_R2[0], MD_TEP_EXCL_B2[0], MD_TEP_EXCL_RC2[0]);
		F.w += MD_TEP_spring_offset[0];
		
		// spring
		number r0 = rmod - MD_TEP_FENE_R0[0];
		number force_module = (MD_TEP_FENE_EPS[0]*r0/(MD_TEP_FENE_DELTA2[0] - SQR(r0)))/rmod;
		F.x += r.x*force_module;
		F.y += r.y*force_module;
		F.z += r.z*force_module;
		F.w += -MD_TEP_FENE_EPS[0]*0.5f*logf(1.f - SQR(r0)/MD_TEP_FENE_DELTA2[0]);

		// bonded twist
		number M = CUDA_DOT(fp, fq) + CUDA_DOT(vp, vq);
		number L = 1.f + CUDA_DOT(up, uq);
		number cos_alpha_plus_gamma = M/L;
		number energy = 0.;
		
		if(-cos_alpha_plus_gamma >= MD_twist_b[0]) printf("We have an issue with the twist angle (%lf) in thread %d\n", cos_alpha_plus_gamma, IND);
		else {
			number4 torque = (MD_kt[0]/L)*(_cross<number, number4>(fp, fq) + _cross<number, number4>(vp, vq) - cos_alpha_plus_gamma*_cross<number, number4>(up, uq))*kt_pref;
			if (-cos_alpha_plus_gamma <= MD_twist_a[0]) {
				energy += MD_kt[0]*(1.f - cos_alpha_plus_gamma);
			}
			else {
				number A = SQR(MD_twist_b[0] - MD_twist_a[0])*CUB(MD_twist_b[0] - MD_twist_a[0])*0.25f;	//A = (b - a)^5/4
				number C = 0.25f*(MD_twist_a[0] - MD_twist_b[0]) + MD_twist_a[0];	// C = (a - b)/4 + a
				torque *= 4.f*A/(SQR(MD_twist_b[0] + cos_alpha_plus_gamma)*CUB(MD_twist_b[0] + cos_alpha_plus_gamma));
				energy += MD_kt[0]*(1.f + A/(CUB(-cos_alpha_plus_gamma - MD_twist_b[0])) + C);
			}
			
			T.x += torque.x;
			T.y += torque.y;
			T.z += torque.z;
		}

		F.w += kt_pref*energy;
		
		// bonded bending
		_bonded_double_bending(up, uq, kb1, kb2, tu, tk, F, T);
		// this is the old "single" bending
		/*number4 torque = MD_kb[0]*kb1*_cross<number, number4>(up, uq);
		T.x += torque.x;
		T.y += torque.y;
		T.z += torque.z;
		F.w += MD_kb[0]*(1.f - CUDA_DOT(up, uq))*kb1;*/
	}
	
	// bonded alignment
	// here we have to take into account the intrinsic directionality of the model (i.e. forces and torques have different expressions if we go from 3' to 5' or vice versa)
	bool keep_r = ((torque_on_p && !alignment_only) || (!torque_on_p && alignment_only));
	number4 ba_tp = (keep_r) ? r : -r;
	number4 ba_up = (torque_on_p) ? up : uq;
	number4 force = MD_ka[0]*(ba_up - ba_tp*CUDA_DOT(ba_up, ba_tp)/SQR(rmod))/rmod;
	if(keep_r) force = -force;
	F.x += force.x;
	F.y += force.y;
	F.z += force.z;
	F.w += MD_ka[0]*(1.f - CUDA_DOT(ba_up, ba_tp)/rmod);
	if(torque_on_p) {
		number4 torque = -(MD_ka[0]*_cross<number, number4>(ba_tp, ba_up))/rmod;
		T.x += torque.x;
		T.y += torque.y;
		T.z += torque.z;
	}
}

template<typename number, typename number4>
__device__ number4 rotate_vector_around_versor(number4 &vector, number4 &versor, number angle){
	number costh = cosf(angle);
	number sinth = sinf(angle);
	number scalar = CUDA_DOT(vector, versor);
	number4 cross = _cross<number, number4>(versor, vector);
	return make_number4<number, number4>(
		versor.x*scalar*(1.f - costh) + vector.x*costh + cross.x*sinth,
		versor.y*scalar*(1.f - costh) + vector.y*costh + cross.y*sinth,
		versor.z*scalar*(1.f - costh) + vector.z*costh + cross.z*sinth,
		0.f
	);
}

template<typename number, typename number4>
__device__ void _twist_boundary_particle(number4 &v1, number4& v3, number4 &T, number4 &o_vect, number4 &w_vect, llint step, number o_modulus) {
	if(MD_twist_boundary_stiff[0] == 0.f) return;
	number4 wt = rotate_vector_around_versor(w_vect, o_vect, step*o_modulus);

	number4 torque_w = MD_twist_boundary_stiff[0]*(_cross<number, number4>(v3, wt));
	number4 torque_o = MD_twist_boundary_stiff[0]*(_cross<number, number4>(v1, o_vect));
	T += torque_w + torque_o;
}

// forces + second step without lists
template <typename number, typename number4>
__global__ void TEP_forces(number4 *poss, GPU_quat<number> *orientations, number4 *forces, number4 *torques,
						LR_bonds *bonds, CUDABox<number, number4> *box, number *kb1_pref, number *kb2_pref,
						number *xk_bending, number *xu_bending, number *kt_pref,
						number4 *o_vects, number4 *w_vects, llint step) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	number4 T = make_number4<number, number4>(0, 0, 0, 0);
	number4 ppos = poss[IND];
	LR_bonds pbonds = bonds[IND];

	number4 a1, a2, a3;
	get_vectors_from_quat<number,number4>(orientations[IND], a1, a2, a3);

	if(pbonds.n3 != P_INVALID) {
		int n3_index = pbonds.n3;
		number4 b1, b2, b3;
		get_vectors_from_quat<number,number4>(orientations[pbonds.n3], b1, b2, b3);

		number kb1 = kb1_pref[n3_index];
		number kb2 = kb2_pref[n3_index];
		number xk = xk_bending[n3_index];
		number xu = xu_bending[n3_index];
		number kt = kt_pref[n3_index];
		_bonded_particle_particle<number, number4, false>(ppos, a1, a2, a3, poss[pbonds.n3], b1, b2, b3, F, T, box, kb1, kb2, xk, xu, kt);
		// if this bead is the last of the strand
		if(pbonds.n5 == P_INVALID) _bonded_particle_particle<number, number4, true>(ppos, a1, a2, a3, poss[pbonds.n3], b1, b2, b3, F, T, box, kb1, kb2, xk, xu, kt, true);
	}
	else _twist_boundary_particle<number, number4>(a1, a3, T, o_vects[0], w_vects[0], step, MD_o_modulus[0]);

	if(pbonds.n5 != P_INVALID) {
		number4 b1, b2, b3;
		get_vectors_from_quat<number,number4>(orientations[pbonds.n5], b1, b2, b3);

		number kb1 = kb1_pref[IND];
		number kb2 = kb2_pref[IND];
		number xk = xk_bending[IND];
		number xu = xu_bending[IND];
		number kt = kt_pref[IND];
		_bonded_particle_particle<number, number4, true>(ppos, a1, a2, a3, poss[pbonds.n5], b1, b2, b3, F, T, box, kb1, kb2, xk, xu, kt);
		// we have to check whether n5 is the last of the strand
		LR_bonds qbonds = bonds[pbonds.n5];
		if(qbonds.n5 == P_INVALID) {
			_bonded_particle_particle<number, number4, false>(ppos, a1, a2, a3, poss[pbonds.n5], b1, b2, b3, F, T, box, kb1, kb2, xk, xu, kt, true);
			_twist_boundary_particle<number, number4>(a1, a3, T, o_vects[1], w_vects[1], step, MD_o_modulus[1]);
		}
	}

	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND && pbonds.n3 != j && pbonds.n5 != j) {
			number4 qpos = poss[j];

			number4 r = box->minimum_image(ppos, qpos);
			_repulsive_lj2(MD_TEP_EXCL_EPS_NONBONDED[0], r, F, MD_TEP_EXCL_S2[0], MD_TEP_EXCL_R2[0], MD_TEP_EXCL_B2[0], MD_TEP_EXCL_RC2[0]);
		}
	}

	F.w *= (number) 0.5f;
	forces[IND] = F;

	T = _vectors_transpose_number4_product(a1, a2, a3, T);
	torques[IND] = T;
}

// forces + second step with verlet lists
template <typename number, typename number4>
__global__ void TEP_forces(number4 *poss, GPU_quat<number> *orientations, number4 *forces, number4 *torques,
						int *matrix_neighs, int *number_neighs, LR_bonds *bonds, CUDABox<number, number4> *box,
						number *kb1_pref, number *kb2_pref,	number *xk_bending, number *xu_bending,
						number *kt_pref, number4 *o_vects, number4 *w_vects, llint step) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	number4 T = make_number4<number, number4>(0, 0, 0, 0);
	number4 ppos = poss[IND];
	LR_bonds pbonds = bonds[IND];

	number4 a1, a2, a3;
	get_vectors_from_quat<number,number4>(orientations[IND], a1, a2, a3);

	if(pbonds.n3 != P_INVALID) {
		int n3_index = pbonds.n3;
		number4 b1, b2, b3;
		get_vectors_from_quat<number,number4>(orientations[pbonds.n3], b1, b2, b3);

		number kb1 = kb1_pref[n3_index];
		number kb2 = kb2_pref[n3_index];
		number xk = xk_bending[n3_index];
		number xu = xu_bending[n3_index];
		number kt = kt_pref[n3_index];
		_bonded_particle_particle<number, number4, false>(ppos, a1, a2, a3, poss[n3_index], b1, b2, b3, F, T, box, kb1, kb2, xk, xu, kt);
		// if this bead is the last of the strand
		if(pbonds.n5 == P_INVALID) _bonded_particle_particle<number, number4, true>(ppos, a1, a2, a3, poss[pbonds.n3], b1, b2, b3, F, T, box, kb1, kb2, xk, xu, kt, true);
	}
	else _twist_boundary_particle<number, number4>(a1, a3, T, o_vects[0], w_vects[0], step, MD_o_modulus[0]);

	if(pbonds.n5 != P_INVALID) {
		number4 b1, b2, b3;
		get_vectors_from_quat<number,number4>(orientations[pbonds.n5], b1, b2, b3);

		number kb1 = kb1_pref[IND];
		number kb2 = kb2_pref[IND];
		number xk = xk_bending[IND];
		number xu = xu_bending[IND];
		number kt = kt_pref[IND];
		_bonded_particle_particle<number, number4, true>(ppos, a1, a2, a3, poss[pbonds.n5], b1, b2, b3, F, T, box, kb1, kb2, xk, xu, kt);
		// we have to check whether n5 is the last of the strand
		LR_bonds qbonds = bonds[pbonds.n5];
		if(qbonds.n5 == P_INVALID) {
			_bonded_particle_particle<number, number4, false>(ppos, a1, a2, a3, poss[pbonds.n5], b1, b2, b3, F, T, box, kb1, kb2, xk, xu, kt, true);
			_twist_boundary_particle<number, number4>(a1, a3, T, o_vects[1], w_vects[1], step, MD_o_modulus[1]);
		}
	}

	int num_neighs = number_neighs[IND];

	for(int j = 0; j < num_neighs; j++) {
		int k_index = matrix_neighs[j*MD_N[0] + IND];

		number4 qpos = poss[k_index];
		number4 r = box->minimum_image(ppos, qpos);
		_repulsive_lj2(MD_TEP_EXCL_EPS_NONBONDED[0], r, F, MD_TEP_EXCL_S2[0], MD_TEP_EXCL_R2[0], MD_TEP_EXCL_B2[0], MD_TEP_EXCL_RC2[0]);
	}

	F.w *= (number) 0.5f;
	forces[IND] = F;

	T = _vectors_transpose_number4_product(a1, a2, a3, T);
	torques[IND] = T;
}

#include "../Lists/CUDASimpleVerletList.h"
#include "../Lists/CUDANoList.h"

template<typename number, typename number4>
CUDATEPInteraction<number, number4>::CUDATEPInteraction() {
	_d_kb1_pref = _d_kb2_pref =_d_kt_pref = _d_xk_bending = _d_xu_bending = NULL;
	_d_o_vects = _d_w_vects = NULL;
	_steps = 0;
}

template<typename number, typename number4>
CUDATEPInteraction<number, number4>::~CUDATEPInteraction() {
	if(_d_kb1_pref != NULL) CUDA_SAFE_CALL( cudaFree(_d_kb1_pref) );
	if(_d_kb2_pref != NULL) CUDA_SAFE_CALL( cudaFree(_d_kb2_pref) );
	if(_d_xk_bending != NULL) CUDA_SAFE_CALL( cudaFree(_d_xk_bending) );
	if(_d_xu_bending != NULL) CUDA_SAFE_CALL( cudaFree(_d_xu_bending) );
	if(_d_kt_pref != NULL) CUDA_SAFE_CALL( cudaFree(_d_kt_pref) );
	if(_d_o_vects != NULL) CUDA_SAFE_CALL( cudaFree(_d_o_vects) );
	if(_d_w_vects != NULL) CUDA_SAFE_CALL( cudaFree(_d_w_vects) );
}

template<typename number, typename number4>
void CUDATEPInteraction<number, number4>::get_settings(input_file &inp) {
	TEPInteraction<number>::get_settings(inp);

	int sort_every;
	if(getInputInt(&inp, "CUDA_sort_every", &sort_every, 0) == KEY_FOUND) {
		if(sort_every > 0) throw oxDNAException("TEP interaction is not compatible with particle sorting, aborting");
	}

	if(this->_prefer_harmonic_over_fene) throw oxDNAException("The 'prefer_harmonic_over_fene' option is not compatible with CUDA");
}

template<typename number, typename number4>
void CUDATEPInteraction<number, number4>::cuda_init(number box_side, int N) {
	CUDABaseInteraction<number, number4>::cuda_init(box_side, N);
	TEPInteraction<number>::init();

	BaseParticle<number> **particles = new BaseParticle<number> *[N];
	int my_N_strands;
	TEPInteraction<number>::read_topology(N, &my_N_strands, particles);

	size_t k_size = N*sizeof(number);
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number>(&_d_kb1_pref, k_size) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number>(&_d_kb2_pref, k_size) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number>(&_d_xk_bending, k_size) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number>(&_d_xu_bending, k_size) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number>(&_d_kt_pref, k_size) );

	CUDA_SAFE_CALL( cudaMemcpy(_d_kb1_pref, this->_kb1_pref, k_size, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_kb2_pref, this->_kb2_pref, k_size, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_xk_bending, this->_xk_bending, k_size, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_xu_bending, this->_xu_bending, k_size, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_kt_pref, this->_kt_pref, k_size, cudaMemcpyHostToDevice) );

	for(int i = 0; i < N; i++) delete particles[i];
	delete[] particles;

	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number4>(&_d_o_vects, 2*sizeof(number4)) );
	CUDA_SAFE_CALL( GpuUtils::LR_cudaMalloc<number4>(&_d_w_vects, 2*sizeof(number4)) );
	number4 h_o_vects[2] = {
			{this->_o1.x, this->_o1.y, this->_o1.z, 0.},
			{this->_o2.x, this->_o2.y, this->_o2.z, 0.}
	};
	number4 h_w_vects[2] = {
			{this->_w1.x, this->_w1.y, this->_w1.z, 0.},
			{this->_w2.x, this->_w2.y, this->_w2.z, 0.}
	};

	CUDA_SAFE_CALL( cudaMemcpy(_d_o_vects, h_o_vects, 2*sizeof(number4), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(_d_w_vects, h_w_vects, 2*sizeof(number4), cudaMemcpyHostToDevice) );
	COPY_NUMBER_TO_FLOAT(MD_twist_boundary_stiff, this->_twist_boundary_stiff);
	float o_modulus[2] = { (float) this->_o1_modulus, (float) this->_o2_modulus };
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_o_modulus, o_modulus, 2*sizeof(float)) );

	CUDA_SAFE_CALL( cudaMemcpyToSymbol(MD_N, &N, sizeof(int)) );

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

template<typename number, typename number4>
void CUDATEPInteraction<number, number4>::compute_forces(CUDABaseList<number, number4> *lists, number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_forces, number4 *d_torques, LR_bonds *d_bonds, CUDABox<number, number4> *d_box) {
	_steps++;

	this->update_increment(_steps);
	this->_time_var = this->update_time_variable(this->_time_var);

	CUDASimpleVerletList<number, number4> *_v_lists = dynamic_cast<CUDASimpleVerletList<number, number4> *>(lists);
	if(_v_lists != NULL) {
		if(_v_lists->use_edge()) throw oxDNAException("use_edge unsupported by TEPInteraction");
		TEP_forces<number, number4>
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_orientations, d_forces, d_torques, _v_lists->_d_matrix_neighs, _v_lists->_d_number_neighs, d_bonds, d_box, _d_kb1_pref, _d_kb2_pref, _d_xk_bending, _d_xu_bending, _d_kt_pref, _d_o_vects, _d_w_vects, this->_time_var);
		CUT_CHECK_ERROR("forces_second_step TEP simple_lists error");
	}

	CUDANoList<number, number4> *_no_lists = dynamic_cast<CUDANoList<number, number4> *>(lists);
	if(_no_lists != NULL) {
		TEP_forces<number, number4>
			<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
			(d_poss, d_orientations, d_forces, d_torques, d_bonds, d_box, _d_kb1_pref, _d_kb2_pref, _d_xk_bending, _d_xu_bending, _d_kt_pref, _d_o_vects, _d_w_vects, this->_time_var);
		CUT_CHECK_ERROR("forces_second_step TEP no_lists error");
	}
}

template class CUDATEPInteraction<float, float4>;
template class CUDATEPInteraction<double, LR_double4>;
