/* System constants */
__constant__ int MD_N[1];
__constant__ int MD_n_forces[1];

__constant__ float MD_hb_multi[1];
__constant__ float MD_F1_A[2];
__constant__ float MD_F1_RC[2];
__constant__ float MD_F1_R0[2];
__constant__ float MD_F1_BLOW[2];
__constant__ float MD_F1_BHIGH[2];
__constant__ float MD_F1_RLOW[2];
__constant__ float MD_F1_RHIGH[2];
__constant__ float MD_F1_RCLOW[2];
__constant__ float MD_F1_RCHIGH[2];
// 50 = 2 * 5 * 5
__constant__ float MD_F1_EPS[50];
__constant__ float MD_F1_SHIFT[50];

__constant__ float MD_F2_K[2];
__constant__ float MD_F2_RC[2];
__constant__ float MD_F2_R0[2];
__constant__ float MD_F2_BLOW[2];
__constant__ float MD_F2_RLOW[2];
__constant__ float MD_F2_RCLOW[2];
__constant__ float MD_F2_BHIGH[2];
__constant__ float MD_F2_RCHIGH[2];
__constant__ float MD_F2_RHIGH[2];

__constant__ float MD_F5_PHI_A[4];
__constant__ float MD_F5_PHI_B[4];
__constant__ float MD_F5_PHI_XC[4];
__constant__ float MD_F5_PHI_XS[4];

__constant__ float MD_dh_RC[1];
__constant__ float MD_dh_RHIGH[1];
__constant__ float MD_dh_prefactor[1];
__constant__ float MD_dh_B[1];
__constant__ float MD_dh_minus_kappa[1];
__constant__ bool MD_dh_half_charged_ends[1];

#include "../cuda_utils/CUDA_lr_common.cuh"

__forceinline__ __device__ void _excluded_volume(const c_number4 &r, c_number4 &F, c_number sigma, c_number rstar, c_number b, c_number rc) {
	c_number rsqr = CUDA_DOT(r, r);

	F.x = F.y = F.z = F.w = (c_number) 0.f;
	if(rsqr < SQR(rc)) {
		if(rsqr > SQR(rstar)) {
			c_number rmod = sqrt(rsqr);
			c_number rrc = rmod - rc;
			c_number fmod = 2.f * EXCL_EPS * b * rrc / rmod;
			F.x = r.x * fmod;
			F.y = r.y * fmod;
			F.z = r.z * fmod;
			F.w = EXCL_EPS * b * SQR(rrc);
		}
		else {
			c_number lj_part = CUB(SQR(sigma)/rsqr);
			c_number fmod = 24.f * EXCL_EPS * (lj_part - 2.f * SQR(lj_part)) / rsqr;
			F.x = r.x * fmod;
			F.y = r.y * fmod;
			F.z = r.z * fmod;
			F.w = 4.f * EXCL_EPS * (SQR(lj_part) - lj_part);
		}
	}
}

__forceinline__ __device__ c_number _f1(c_number r, int type, int n3, int n5) {
	c_number val = (c_number) 0.f;
	if(r < MD_F1_RCHIGH[type]) {
		int eps_index = 25 * type + n3 * 5 + n5;
		if(r > MD_F1_RHIGH[type]) {
			val = MD_F1_EPS[eps_index] * MD_F1_BHIGH[type] * SQR(r - MD_F1_RCHIGH[type]);
		}
		else if(r > MD_F1_RLOW[type]) {
			c_number tmp = 1.f - expf(-(r - MD_F1_R0[type]) * MD_F1_A[type]);
			val = MD_F1_EPS[eps_index] * SQR(tmp) - MD_F1_SHIFT[eps_index];
		}
		else if(r > MD_F1_RCLOW[type]) {
			val = MD_F1_EPS[eps_index] * MD_F1_BLOW[type] * SQR(r - MD_F1_RCLOW[type]);
		}
	}

	return val;
}

__forceinline__ __device__ c_number _f1D(c_number r, int type, int n3, int n5) {
	c_number val = (c_number) 0.f;
	if(r < MD_F1_RCHIGH[type]) {
		float eps = MD_F1_EPS[25 * type + n3 * 5 + n5];
		if(r > MD_F1_RHIGH[type]) {
			val = 2.f * eps * MD_F1_BHIGH[type] * (r - MD_F1_RCHIGH[type]);
		}
		else if(r > MD_F1_RLOW[type]) {
			c_number tmp = expf(-(r - MD_F1_R0[type]) * MD_F1_A[type]);
			val = 2.f * eps * (1.f - tmp) * tmp * MD_F1_A[type];
		}
		else if(r > MD_F1_RCLOW[type]) {
			val = 2.f * eps * MD_F1_BLOW[type] * (r - MD_F1_RCLOW[type]);
		}
	}

	return val;
}

__forceinline__ __device__ c_number _f2(c_number r, int type) {
	c_number val = (c_number) 0.f;
	if(r < MD_F2_RCHIGH[type]) {
		if(r > MD_F2_RHIGH[type]) {
			val = MD_F2_K[type] * MD_F2_BHIGH[type] * SQR(r - MD_F2_RCHIGH[type]);
		}
		else if(r > MD_F2_RLOW[type]) {
			val = (MD_F2_K[type] * 0.5f) * (SQR(r - MD_F2_R0[type]) - SQR(MD_F2_RC[type] - MD_F2_R0[type]));
		}
		else if(r > MD_F2_RCLOW[type]) {
			val = MD_F2_K[type] * MD_F2_BLOW[type] * SQR(r - MD_F2_RCLOW[type]);
		}
	}
	return val;
}

__forceinline__ __device__ c_number _f2D(c_number r, int type) {
	c_number val = (c_number) 0.f;
	if(r < MD_F2_RCHIGH[type]) {
		if(r > MD_F2_RHIGH[type]) {
			val = 2.f * MD_F2_K[type] * MD_F2_BHIGH[type] * (r - MD_F2_RCHIGH[type]);
		}
		else if(r > MD_F2_RLOW[type]) {
			val = MD_F2_K[type] * (r - MD_F2_R0[type]);
		}
		else if(r > MD_F2_RCLOW[type]) {
			val = 2.f * MD_F2_K[type] * MD_F2_BLOW[type] * (r - MD_F2_RCLOW[type]);
		}
	}
	return val;
}

__forceinline__ __device__ c_number _f4(c_number t, float t0, float ts, float tc, float a, float b) {
	c_number val = (c_number) 0.f;
	t = copysignf(t - t0, (c_number) 1.f);

	if(t < tc) {
		val = (t > ts) ? b * SQR(tc - t) : 1.f - a * SQR(t);
	}

	return val;
}

__forceinline__ __device__ c_number _f4_pure_harmonic(c_number t, float a, float b) {
	// for getting a f4t1 function with a continuous derivative that is less disruptive to the potential
	t -= b;
	return (t < 0) ? 0.f : a * SQR(t);
}

__forceinline__ __device__ c_number _f4D(c_number t, float t0, float ts, float tc, float a, float b) {
	c_number val = (c_number) 0.f;
	t -= t0;
	// this function is a parabola centered in t0. If tt0 < 0 then the value of the function
	// is the same but the value of its derivative has the opposite sign, so m = -1
	c_number m = copysignf((c_number) 1.f, t);
	t = copysignf(t, (c_number) 1.f);

	if(t < tc) {
		val = (t > ts) ? 2.f * m * b * (t - tc) : -2.f * m * a * t;
	}

	return val;
}

__forceinline__ __device__ c_number _f4D_pure_harmonic(c_number t, float a, float b) {
	// for getting a f4t1 function with a continuous derivative that is less disruptive to the potential
	t -= b;
	return (t < 0) ? 0.f : 2.f * a * t;
}

__forceinline__ __device__ c_number _f5(c_number f, int type) {
	c_number val = (c_number) 0.f;

	if(f > MD_F5_PHI_XC[type]) {
		if(f < MD_F5_PHI_XS[type]) {
			val = MD_F5_PHI_B[type] * SQR(MD_F5_PHI_XC[type] - f);
		}
		else if(f < 0.f) {
			val = (c_number) 1.f - MD_F5_PHI_A[type] * SQR(f);
		}
		else val = 1.f;
	}

	return val;
}

__forceinline__ __device__ c_number _f5D(c_number f, int type) {
	c_number val = (c_number) 0.f;

	if(f > MD_F5_PHI_XC[type]) {
		if(f < MD_F5_PHI_XS[type]) {
			val = 2.f * MD_F5_PHI_B[type] * (f - MD_F5_PHI_XC[type]);
		}
		else if(f < 0.f) {
			val = (c_number) -2.f * MD_F5_PHI_A[type] * f;
		}
	}

	return val;
}

template<bool qIsN3>
__device__ void _bonded_excluded_volume(const c_number4 &r, const c_number4 &n3pos_base, const c_number4 &n3pos_back, const c_number4 &n5pos_base,
		const c_number4 &n5pos_back, c_number4 &F, c_number4 &T) {
	c_number4 Ftmp;
	// BASE-BASE
	c_number4 rcenter = r + n3pos_base - n5pos_base;
	_excluded_volume(rcenter, Ftmp, EXCL_S2, EXCL_R2, EXCL_B2, EXCL_RC2);
	c_number4 torquep1 = (qIsN3) ? _cross(n5pos_base, Ftmp) : _cross(n3pos_base, Ftmp);
	F += Ftmp;

	// n5-BASE vs. n3-BACK
	rcenter = r + n3pos_back - n5pos_base;
	_excluded_volume(rcenter, Ftmp, EXCL_S3, EXCL_R3, EXCL_B3, EXCL_RC3);
	c_number4 torquep2 = (qIsN3) ? _cross(n5pos_base, Ftmp) : _cross(n3pos_back, Ftmp);
	F += Ftmp;

	// n5-BACK vs. n3-BASE
	rcenter = r + n3pos_base - n5pos_back;
	_excluded_volume(rcenter, Ftmp, EXCL_S4, EXCL_R4, EXCL_B4, EXCL_RC4);
	c_number4 torquep3 = (qIsN3) ? _cross(n5pos_back, Ftmp) : _cross(n3pos_base, Ftmp);
	F += Ftmp;

	T += torquep1 + torquep2 + torquep3;
}

template<bool qIsN3>
__device__ void _bonded_part(const c_number4 &r, const c_number4 &n5pos, const c_number4 &n5x, const c_number4 &n5y, const c_number4 &n5z, const c_number4 &n3pos, const c_number4 &n3x,
		const c_number4 &n3y, const c_number4 &n3z, c_number4 &F, c_number4 &T, bool grooving, bool use_oxDNA2_FENE, bool use_mbf,
		c_number mbf_xmax, c_number mbf_finf) {

	int n3type = get_particle_type(n3pos);
	int n5type = get_particle_type(n5pos);

	c_number4 n5pos_back;
	if(grooving) n5pos_back = n5x * POS_MM_BACK1 + n5y * POS_MM_BACK2;
	else n5pos_back = n5x * POS_BACK;

	c_number4 n5pos_base = n5x * POS_BASE;
	c_number4 n5pos_stack = n5x * POS_STACK;

	c_number4 n3pos_back;
	if(grooving) n3pos_back = n3x * POS_MM_BACK1 + n3y * POS_MM_BACK2;
	else n3pos_back = n3x * POS_BACK;
	c_number4 n3pos_base = n3x * POS_BASE;
	c_number4 n3pos_stack = n3x * POS_STACK;

	c_number4 rback = r + n3pos_back - n5pos_back;
	c_number rbackmod = _module(rback);
	c_number rbackr0;
	if(use_oxDNA2_FENE) rbackr0 = rbackmod - FENE_R0_OXDNA2;
	else rbackr0 = rbackmod - FENE_R0_OXDNA;

	c_number4 Ftmp;
	if(use_mbf == true && fabsf(rbackr0) > mbf_xmax) {
		// this is the "relax" potential, i.e. the standard FENE up to xmax and then something like A + B log(r) for r>xmax
		c_number fene_xmax = -(FENE_EPS / 2.f) * logf(1.f - mbf_xmax * mbf_xmax / FENE_DELTA2);
		c_number mbf_fmax = (FENE_EPS * mbf_xmax / (FENE_DELTA2 - SQR(mbf_xmax)));
		c_number long_xmax = (mbf_fmax - mbf_finf) * mbf_xmax * logf(mbf_xmax) + mbf_finf * mbf_xmax;
		Ftmp = rback * (copysignf(1.f, rbackr0) * ((mbf_fmax - mbf_finf) * mbf_xmax / fabsf(rbackr0) + mbf_finf) / rbackmod);
		Ftmp.w = (mbf_fmax - mbf_finf) * mbf_xmax * logf(fabsf(rbackr0)) + mbf_finf * fabsf(rbackr0) - long_xmax + fene_xmax;
	}
	else {
		Ftmp = rback * ((FENE_EPS * rbackr0 / (FENE_DELTA2 - SQR(rbackr0))) / rbackmod);
		Ftmp.w = -FENE_EPS * ((c_number) 0.5f) * logf(1 - SQR(rbackr0) / FENE_DELTA2);
	}

	c_number4 Ttmp = (qIsN3) ? _cross(n5pos_back, Ftmp) : _cross(n3pos_back, Ftmp);
	// EXCLUDED VOLUME
	_bonded_excluded_volume<qIsN3>(r, n3pos_base, n3pos_back, n5pos_base, n5pos_back, Ftmp, Ttmp);

	if(qIsN3) {
		F += Ftmp;
		T += Ttmp;
	}
	else {
		F -= Ftmp;
		T -= Ttmp;
	}

	// STACKING
	c_number4 rstack = r + n3pos_stack - n5pos_stack;
	c_number rstackmod = _module(rstack);
	c_number4 rstackdir = make_c_number4(rstack.x / rstackmod, rstack.y / rstackmod, rstack.z / rstackmod, 0);
	// This is the position the backbone would have with major-minor grooves the same width.
	// We need to do this to implement different major-minor groove widths because rback is
	// used as a reference point for things that have nothing to do with the actual backbone
	// position (in this case, the stacking interaction).
	c_number4 rbackref = r + n3x * POS_BACK - n5x * POS_BACK;
	c_number rbackrefmod = _module(rbackref);

	c_number t4 = CUDA_LRACOS(CUDA_DOT(n3z, n5z));
	c_number cost5 = CUDA_DOT(n5z, rstackdir);
	c_number t5 = CUDA_LRACOS(cost5);
	c_number cost6 = -CUDA_DOT(n3z, rstackdir);
	c_number t6 = CUDA_LRACOS(cost6);
	c_number cosphi1 = CUDA_DOT(n5y, rbackref) / rbackrefmod;
	c_number cosphi2 = CUDA_DOT(n3y, rbackref) / rbackrefmod;

	// functions
	c_number f1 = _f1(rstackmod, STCK_F1, n3type, n5type);
	c_number f4t4 = _f4(t4, STCK_THETA4_T0, STCK_THETA4_TS, STCK_THETA4_TC, STCK_THETA4_A, STCK_THETA4_B);
	c_number f4t5 = _f4(PI - t5, STCK_THETA5_T0, STCK_THETA5_TS, STCK_THETA5_TC, STCK_THETA5_A, STCK_THETA5_B);
	c_number f4t6 = _f4(t6, STCK_THETA6_T0, STCK_THETA6_TS, STCK_THETA6_TC, STCK_THETA6_A, STCK_THETA6_B);
	c_number f5phi1 = _f5(cosphi1, STCK_F5_PHI1);
	c_number f5phi2 = _f5(cosphi2, STCK_F5_PHI2);

	c_number energy = f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;

	if(energy != (c_number) 0) {
		// and their derivatives
		c_number f1D = _f1D(rstackmod, STCK_F1, n3type, n5type);
		c_number f4t4D = _f4D(t4, STCK_THETA4_T0, STCK_THETA4_TS, STCK_THETA4_TC, STCK_THETA4_A, STCK_THETA4_B);
		c_number f4t5D = _f4D(PI - t5, STCK_THETA5_T0, STCK_THETA5_TS, STCK_THETA5_TC, STCK_THETA5_A, STCK_THETA5_B);
		c_number f4t6D = _f4D(t6, STCK_THETA6_T0, STCK_THETA6_TS, STCK_THETA6_TC, STCK_THETA6_A, STCK_THETA6_B);
		c_number f5phi1D = _f5D(cosphi1, STCK_F5_PHI1);
		c_number f5phi2D = _f5D(cosphi2, STCK_F5_PHI2);

		// RADIAL
		Ftmp = rstackdir * (energy * f1D / f1);

		// THETA 5
		Ftmp += stably_normalised(n5z - cost5 * rstackdir) * (energy * f4t5D / (f4t5 * rstackmod));

		// THETA 6
		Ftmp += stably_normalised(n3z + cost6 * rstackdir) * (energy * f4t6D / (f4t6 * rstackmod));

		// COS PHI 1
		// here particle p is referred to using the a while particle q is referred with the b
		c_number ra2 = CUDA_DOT(rstackdir, n5y);
		c_number ra1 = CUDA_DOT(rstackdir, n5x);
		c_number rb1 = CUDA_DOT(rstackdir, n3x);
		c_number a2b1 = CUDA_DOT(n5y, n3x);
		c_number dcosphi1dr = (SQR(rstackmod) * ra2 - ra2 * SQR(rbackrefmod) - rstackmod * (a2b1 + ra2 * (-ra1 + rb1)) * GAMMA + a2b1 * (-ra1 + rb1) * SQR(GAMMA)) / (SQR(rbackrefmod) * rbackrefmod);
		c_number dcosphi1dra1 = rstackmod * GAMMA * (rstackmod * ra2 - a2b1 * GAMMA) / (SQR(rbackrefmod) * rbackrefmod);
		c_number dcosphi1dra2 = -rstackmod / rbackrefmod;
		c_number dcosphi1drb1 = -(rstackmod * GAMMA * (rstackmod * ra2 - a2b1 * GAMMA)) / (SQR(rbackrefmod) * rbackrefmod);
		c_number dcosphi1da1b1 = SQR(GAMMA) * (-rstackmod * ra2 + a2b1 * GAMMA) / (SQR(rbackrefmod) * rbackrefmod);
		c_number dcosphi1da2b1 = GAMMA / rbackrefmod;

		c_number force_part_phi1 = energy * f5phi1D / f5phi1;

		Ftmp -= (rstackdir * dcosphi1dr + ((n5y - ra2 * rstackdir) * dcosphi1dra2 + (n5x - ra1 * rstackdir) * dcosphi1dra1 + (n3x - rb1 * rstackdir) * dcosphi1drb1) / rstackmod) * force_part_phi1;

		// COS PHI 2
		// here particle p -> b, particle q -> a
		ra2 = CUDA_DOT(rstackdir, n3y);
		ra1 = rb1;
		rb1 = CUDA_DOT(rstackdir, n5x);
		a2b1 = CUDA_DOT(n3y, n5x);
		c_number dcosphi2dr = ((rstackmod * ra2 + a2b1 * GAMMA) * (rstackmod + (rb1 - ra1) * GAMMA) - ra2 * SQR(rbackrefmod)) / (SQR(rbackrefmod) * rbackrefmod);
		c_number dcosphi2dra1 = -rstackmod * GAMMA * (rstackmod * ra2 + a2b1 * GAMMA) / (SQR(rbackrefmod) * rbackrefmod);
		c_number dcosphi2dra2 = -rstackmod / rbackrefmod;
		c_number dcosphi2drb1 = (rstackmod * GAMMA * (rstackmod * ra2 + a2b1 * GAMMA)) / (SQR(rbackrefmod) * rbackrefmod);
		c_number dcosphi2da1b1 = -SQR(GAMMA) * (rstackmod * ra2 + a2b1 * GAMMA) / (SQR(rbackrefmod) * rbackrefmod);
		c_number dcosphi2da2b1 = -GAMMA / rbackrefmod;

		c_number force_part_phi2 = energy * f5phi2D / f5phi2;

		Ftmp -= (rstackdir * dcosphi2dr + ((n3y - rstackdir * ra2) * dcosphi2dra2 + (n3x - rstackdir * ra1) * dcosphi2dra1 + (n5x - rstackdir * rb1) * dcosphi2drb1) / rstackmod) * force_part_phi2;

		if(qIsN3) Ttmp = _cross(n5pos_stack, Ftmp);
		else Ttmp = _cross(n3pos_stack, Ftmp);

		// THETA 4
		Ttmp += stably_normalised(_cross(n3z, n5z)) * (-energy * f4t4D / f4t4);

		// PHI 1 & PHI 2
		if(qIsN3) {
			Ttmp += (-force_part_phi1 * dcosphi1dra2) * _cross(rstackdir, n5y) - _cross(rstackdir, n5x) * force_part_phi1 * dcosphi1dra1;

			Ttmp += (-force_part_phi2 * dcosphi2drb1) * _cross(rstackdir, n5x);
		}
		else {
			Ttmp += force_part_phi1 * dcosphi1drb1 * _cross(rstackdir, n3x);

			Ttmp += force_part_phi2 * dcosphi2dra2 * _cross(rstackdir, n3y) + force_part_phi2 * dcosphi2dra1 * _cross(rstackdir, n3x);
		}

		Ttmp += force_part_phi1 * dcosphi1da2b1 * _cross(n5y, n3x) + _cross(n5x, n3x) * force_part_phi1 * dcosphi1da1b1;

		Ttmp += force_part_phi2 * dcosphi2da2b1 * _cross(n5x, n3y) + _cross(n5x, n3x) * force_part_phi2 * dcosphi2da1b1;

		Ftmp.w = energy;
		if(qIsN3) {
			// THETA 5
			Ttmp += stably_normalised(_cross(rstackdir, n5z)) * energy * f4t5D / f4t5;

			T += Ttmp;
			F += Ftmp;
		}
		else {
			// THETA 6
			Ttmp += stably_normalised(_cross(rstackdir, n3z)) * (-energy * f4t6D / f4t6);

			T -= Ttmp;
			F -= Ftmp;
		}
	}
}

__device__ void _particle_particle_DNA_interaction(const c_number4 &r, const c_number4 &ppos, const c_number4 &a1, const c_number4 &a2, const c_number4 &a3,
		const c_number4 &qpos, const c_number4 &b1,	const c_number4 &b2, const c_number4 &b3, c_number4 &F, c_number4 &T, bool grooving,
		bool use_debye_huckel, bool use_oxDNA2_coaxial_stacking,
		bool p_is_end, bool q_is_end) {
	int ptype = get_particle_type(ppos);
	int qtype = get_particle_type(qpos);
	int pbtype = get_particle_btype(ppos);
	int qbtype = get_particle_btype(qpos);
	int int_type = pbtype + qbtype;

	c_number4 ppos_back = (grooving) ? POS_MM_BACK1 * a1 + POS_MM_BACK2 * a2 : POS_BACK * a1;
	c_number4 ppos_base = POS_BASE * a1;
	c_number4 ppos_stack = POS_STACK * a1;

	c_number4 qpos_back = (grooving) ? POS_MM_BACK1 * b1 + POS_MM_BACK2 * b2 : POS_BACK * b1;
	c_number4 qpos_base = POS_BASE * b1;
	c_number4 qpos_stack = POS_STACK * b1;

	c_number tot_Tw = T.w;

	// excluded volume
	// BACK-BACK
	c_number4 Ftmp = make_c_number4(0, 0, 0, 0);
	c_number4 rbackbone = r + qpos_back - ppos_back;
	_excluded_volume(rbackbone, Ftmp, EXCL_S1, EXCL_R1, EXCL_B1, EXCL_RC1);
	c_number4 Ttmp = _cross(ppos_back, Ftmp);
	_bonded_excluded_volume<true>(r, qpos_base, qpos_back, ppos_base, ppos_back, Ftmp, Ttmp);

	F += Ftmp;

	// DEBYE HUCKEL
	if(use_debye_huckel) {
		c_number rbackmod = _module(rbackbone);
		if(rbackmod < MD_dh_RC[0]) {
			c_number4 rbackdir = rbackbone / rbackmod;
			if(rbackmod < MD_dh_RHIGH[0]) {
				Ftmp = rbackdir * (-MD_dh_prefactor[0] * expf(MD_dh_minus_kappa[0] * rbackmod) * (MD_dh_minus_kappa[0] / rbackmod - 1.0f / SQR(rbackmod)));
				Ftmp.w = expf(rbackmod * MD_dh_minus_kappa[0]) * (MD_dh_prefactor[0] / rbackmod);
			}
			else {
				Ftmp = rbackdir * (-2.0f * MD_dh_B[0] * (rbackmod - MD_dh_RC[0]));
				Ftmp.w = MD_dh_B[0] * SQR(rbackmod - MD_dh_RC[0]);
			}

			// check for half-charge strand ends
			if(MD_dh_half_charged_ends[0] && p_is_end) {
				Ftmp *= 0.5f;
			}
			if(MD_dh_half_charged_ends[0] && q_is_end) {
				Ftmp *= 0.5f;
			}

			Ttmp -= _cross(ppos_back, Ftmp);
			F -= Ftmp;
		}
	}

	// HYDROGEN BONDING
	c_number4 rhydro = r + qpos_base - ppos_base;
	c_number rhydromodsqr = CUDA_DOT(rhydro, rhydro);
	if(int_type == 3 && SQR(HYDR_RCLOW) < rhydromodsqr && rhydromodsqr < SQR(HYDR_RCHIGH)) {
		c_number hb_multi = (abs(qbtype) >= 300 && abs(pbtype) >= 300) ? MD_hb_multi[0] : 1.f;
		// versor and magnitude of the base-base separation
		c_number rhydromod = sqrtf(rhydromodsqr);
		c_number4 rhydrodir = rhydro / rhydromod;

		// angles involved in the HB interaction
		c_number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));
		c_number cost2 = -CUDA_DOT(b1, rhydrodir);
		c_number t2 = CUDA_LRACOS(cost2);
		c_number cost3 = CUDA_DOT(a1, rhydrodir);
		c_number t3 = CUDA_LRACOS(cost3);
		c_number t4 = CUDA_LRACOS(CUDA_DOT(a3, b3));
		c_number cost7 = -CUDA_DOT(rhydrodir, b3);
		c_number t7 = CUDA_LRACOS(cost7);
		c_number cost8 = CUDA_DOT(rhydrodir, a3);
		c_number t8 = CUDA_LRACOS(cost8);

		// functions called at their relevant arguments
		c_number f1 = hb_multi * _f1(rhydromod, HYDR_F1, ptype, qtype);
		c_number f4t1 = _f4(t1, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC, HYDR_THETA1_A, HYDR_THETA1_B);
		c_number f4t2 = _f4(t2, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC, HYDR_THETA2_A, HYDR_THETA2_B);
		c_number f4t3 = _f4(t3, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC, HYDR_THETA3_A, HYDR_THETA3_B);
		c_number f4t4 = _f4(t4, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC, HYDR_THETA4_A, HYDR_THETA4_B);
		c_number f4t7 = _f4(t7, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC, HYDR_THETA7_A, HYDR_THETA7_B);
		c_number f4t8 = _f4(t8, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC, HYDR_THETA8_A, HYDR_THETA8_B);

		c_number hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
		tot_Tw += hb_energy;

		if(hb_energy < (c_number) 0) {
			// derivatives called at the relevant arguments
			c_number f1D = hb_multi * _f1D(rhydromod, HYDR_F1, ptype, qtype);
			c_number f4t1D = -_f4D(t1, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC, HYDR_THETA1_A, HYDR_THETA1_B);
			c_number f4t2D = -_f4D(t2, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC, HYDR_THETA2_A, HYDR_THETA2_B);
			c_number f4t3D = _f4D(t3, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC, HYDR_THETA3_A, HYDR_THETA3_B);
			c_number f4t4D = _f4D(t4, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC, HYDR_THETA4_A, HYDR_THETA4_B);
			c_number f4t7D = -_f4D(t7, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC, HYDR_THETA7_A, HYDR_THETA7_B);
			c_number f4t8D = _f4D(t8, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC, HYDR_THETA8_A, HYDR_THETA8_B);

			// RADIAL PART
			Ftmp = rhydrodir * hb_energy * f1D / f1;

			// TETA4; t4 = LRACOS (a3 * b3);
			Ttmp -= stably_normalised(_cross(a3, b3)) * (-hb_energy * f4t4D / f4t4);

			// TETA1; t1 = LRACOS (-a1 * b1);
			Ttmp -= stably_normalised(_cross(a1, b1)) * (-hb_energy * f4t1D / f4t1);

			// TETA2; t2 = LRACOS (-b1 * rhydrodir);
			Ftmp -= stably_normalised(b1 + rhydrodir * cost2) * (hb_energy * f4t2D / (f4t2 * rhydromod));

			// TETA3; t3 = LRACOS (a1 * rhydrodir);
			c_number part = -hb_energy * f4t3D / f4t3;
			Ftmp -= stably_normalised(a1 - rhydrodir * cost3) * (-part / rhydromod);
			Ttmp += stably_normalised(_cross(rhydrodir, a1)) * part;

			// THETA7; t7 = LRACOS (-rhydrodir * b3);
			Ftmp -= stably_normalised(b3 + rhydrodir * cost7) * (hb_energy * f4t7D / (f4t7 * rhydromod));

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			part = -hb_energy * f4t8D / f4t8;
			Ftmp -= stably_normalised(a3 - rhydrodir * cost8) * (-part / rhydromod);
			Ttmp += stably_normalised(_cross(rhydrodir, a3)) * part;

			Ttmp += _cross(ppos_base, Ftmp);

			Ftmp.w = hb_energy;
			F += Ftmp;
		}
	}
	// END HYDROGEN BONDING

	// CROSS STACKING
	c_number4 rcstack = rhydro;
	c_number rcstackmodsqr = rhydromodsqr;
	if(SQR(CRST_RCLOW) < rcstackmodsqr && rcstackmodsqr < SQR(CRST_RCHIGH)) {
		c_number rcstackmod = sqrtf(rcstackmodsqr);
		c_number4 rcstackdir = rcstack / rcstackmod;

		// angles involved in the CSTCK interaction
		c_number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));
		c_number cost2 = -CUDA_DOT(b1, rcstackdir);
		c_number t2 = CUDA_LRACOS(cost2);
		c_number cost3 = CUDA_DOT(a1, rcstackdir);
		c_number t3 = CUDA_LRACOS(cost3);
		c_number t4 = CUDA_LRACOS(CUDA_DOT(a3, b3));
		c_number cost7 = -CUDA_DOT(rcstackdir, b3);
		c_number t7 = CUDA_LRACOS(cost7);
		c_number cost8 = CUDA_DOT(rcstackdir, a3);
		c_number t8 = CUDA_LRACOS(cost8);

		// functions called at their relevant arguments
		c_number f2 = _f2(rcstackmod, CRST_F2);
		c_number f4t1 = _f4(t1, CRST_THETA1_T0, CRST_THETA1_TS, CRST_THETA1_TC, CRST_THETA1_A, CRST_THETA1_B);
		c_number f4t2 = _f4(t2, CRST_THETA2_T0, CRST_THETA2_TS, CRST_THETA2_TC, CRST_THETA2_A, CRST_THETA2_B);
		c_number f4t3 = _f4(t3, CRST_THETA3_T0, CRST_THETA3_TS, CRST_THETA3_TC, CRST_THETA3_A, CRST_THETA3_B);
		c_number f4t4 = _f4(t4, CRST_THETA4_T0, CRST_THETA4_TS, CRST_THETA4_TC, CRST_THETA4_A, CRST_THETA4_B) + _f4(PI - t4, CRST_THETA4_T0, CRST_THETA4_TS, CRST_THETA4_TC, CRST_THETA4_A, CRST_THETA4_B);
		c_number f4t7 = _f4(t7, CRST_THETA7_T0, CRST_THETA7_TS, CRST_THETA7_TC, CRST_THETA7_A, CRST_THETA7_B) + _f4(PI - t7, CRST_THETA7_T0, CRST_THETA7_TS, CRST_THETA7_TC, CRST_THETA7_A, CRST_THETA7_B);
		c_number f4t8 = _f4(t8, CRST_THETA8_T0, CRST_THETA8_TS, CRST_THETA8_TC, CRST_THETA8_A, CRST_THETA8_B) + _f4(PI - t8, CRST_THETA8_T0, CRST_THETA8_TS, CRST_THETA8_TC, CRST_THETA8_A, CRST_THETA8_B);

		c_number cstk_energy = f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

		if(cstk_energy < (c_number) 0) {
			// derivatives called at the relevant arguments
			c_number f2D = _f2D(rcstackmod, CRST_F2);
			c_number f4t1D = -_f4D(t1, CRST_THETA1_T0, CRST_THETA1_TS, CRST_THETA1_TC, CRST_THETA1_A, CRST_THETA1_B);
			c_number f4t2D = -_f4D(t2, CRST_THETA2_T0, CRST_THETA2_TS, CRST_THETA2_TC, CRST_THETA2_A, CRST_THETA2_B);
			c_number f4t3D = _f4D(t3, CRST_THETA3_T0, CRST_THETA3_TS, CRST_THETA3_TC, CRST_THETA3_A, CRST_THETA3_B);
			c_number f4t4D = _f4D(t4, CRST_THETA4_T0, CRST_THETA4_TS, CRST_THETA4_TC, CRST_THETA4_A, CRST_THETA4_B) - _f4D(PI - t4, CRST_THETA4_T0, CRST_THETA4_TS, CRST_THETA4_TC, CRST_THETA4_A, CRST_THETA4_B);
			c_number f4t7D = -_f4D(t7, CRST_THETA7_T0, CRST_THETA7_TS, CRST_THETA7_TC, CRST_THETA7_A, CRST_THETA7_B) + _f4D(PI - t7, CRST_THETA7_T0, CRST_THETA7_TS, CRST_THETA7_TC, CRST_THETA7_A, CRST_THETA7_B);
			c_number f4t8D = _f4D(t8, CRST_THETA8_T0, CRST_THETA8_TS, CRST_THETA8_TC, CRST_THETA8_A, CRST_THETA8_B) - _f4D(PI - t8, CRST_THETA8_T0, CRST_THETA8_TS, CRST_THETA8_TC, CRST_THETA8_A, CRST_THETA8_B);

			// RADIAL PART
			Ftmp = rcstackdir * (cstk_energy * f2D / f2);

			// THETA1; t1 = LRACOS (-a1 * b1);
			Ttmp -= stably_normalised(_cross(a1, b1)) * (-cstk_energy * f4t1D / f4t1);

			// TETA2; t2 = LRACOS (-b1 * rhydrodir);
			Ftmp -= stably_normalised(b1 + rcstackdir * cost2) * (cstk_energy * f4t2D / (f4t2 * rcstackmod));

			// TETA3; t3 = LRACOS (a1 * rhydrodir);
			c_number part = -cstk_energy * f4t3D / f4t3;
			Ftmp -= stably_normalised(a1 - rcstackdir * cost3) * (-part / rcstackmod);
			Ttmp += stably_normalised(_cross(rcstackdir, a1)) * part;

			// TETA4; t4 = LRACOS (a3 * b3);
			Ttmp -= stably_normalised(_cross(a3, b3)) * (-cstk_energy * f4t4D / f4t4);

			// THETA7; t7 = LRACOS (-rcsrackir * b3);
			Ftmp -= stably_normalised(b3 + rcstackdir * cost7) * (cstk_energy * f4t7D / (f4t7 * rcstackmod));

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			part = -cstk_energy * f4t8D / f4t8;
			Ftmp -= stably_normalised(a3 - rcstackdir * cost8) * (-part / rcstackmod);
			Ttmp += stably_normalised(_cross(rcstackdir, a3)) * part;

			Ttmp += _cross(ppos_base, Ftmp);

			Ftmp.w = cstk_energy;
			F += Ftmp;
		}
	}

	// COAXIAL STACKING
	// oxDNA1's coaxial stacking was a bit different in terms of parameters, and also had an additional term
	// here we try to avoid duplications by using if statements (which are always true or false and therefore
	// do not branch out), or use ternary operators (sparsely, since they decrease the readability)
	c_number4 rstack = r + qpos_stack - ppos_stack;
	c_number rstackmodsqr = CUDA_DOT(rstack, rstack);
	if(SQR(CXST_RCLOW) < rstackmodsqr && rstackmodsqr < SQR(CXST_RCHIGH)) {
		c_number rstackmod = sqrtf(rstackmodsqr);
		c_number4 rstackdir = rstack / rstackmod;

		// angles involved in the CXST interaction
		c_number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));
		c_number t4 = CUDA_LRACOS(CUDA_DOT(a3, b3));
		c_number cost5 = CUDA_DOT(a3, rstackdir);
		c_number t5 = CUDA_LRACOS(cost5);
		c_number cost6 = -CUDA_DOT(b3, rstackdir);
		c_number t6 = CUDA_LRACOS(cost6);

		c_number cosphi3 = 1.f;
		c_number f5cosphi3 = 1.f;
		c_number rbackrefmod = 1.f;
		if(!use_oxDNA2_coaxial_stacking) {
			// This is the position the backbone would have with major-minor grooves the same width.
			// We need to do this to implement different major-minor groove widths because rback is
			// used as a reference point for things that have nothing to do with the actual backbone
			// position (in this case, the coaxial stacking interaction).
			c_number4 rbackboneref = r + POS_BACK * b1 - POS_BACK * a1;
			rbackrefmod = _module(rbackboneref);
			c_number4 rbackbonerefdir = rbackboneref / rbackrefmod;
			cosphi3 = CUDA_DOT(rstackdir, (_cross(rbackbonerefdir, a1)));
			f5cosphi3 = _f5(cosphi3, CXST_F5_PHI3);
		}

		// functions called at their relevant arguments
		c_number f2 = _f2(rstackmod, CXST_F2);
		c_number f4t1 = (use_oxDNA2_coaxial_stacking) ?
				_f4(t1, CXST_THETA1_T0_OXDNA2, CXST_THETA1_TS, CXST_THETA1_TC, CXST_THETA1_A, CXST_THETA1_B) + _f4_pure_harmonic(t1, CXST_THETA1_SA, CXST_THETA1_SB) :
				_f4(t1, CXST_THETA1_T0_OXDNA, CXST_THETA1_TS, CXST_THETA1_TC, CXST_THETA1_A, CXST_THETA1_B) + _f4(2 * PI - t1, CXST_THETA1_T0_OXDNA, CXST_THETA1_TS, CXST_THETA1_TC, CXST_THETA1_A, CXST_THETA1_B);
		c_number f4t4 = _f4(t4, CXST_THETA4_T0, CXST_THETA4_TS, CXST_THETA4_TC, CXST_THETA4_A, CXST_THETA4_B);
		c_number f4t5 = _f4(t5, CXST_THETA5_T0, CXST_THETA5_TS, CXST_THETA5_TC, CXST_THETA5_A, CXST_THETA5_B) + _f4(PI - t5, CXST_THETA5_T0, CXST_THETA5_TS, CXST_THETA5_TC, CXST_THETA5_A, CXST_THETA5_B);
		c_number f4t6 = _f4(t6, CXST_THETA6_T0, CXST_THETA6_TS, CXST_THETA6_TC, CXST_THETA6_A, CXST_THETA6_B) + _f4(PI - t6, CXST_THETA6_T0, CXST_THETA6_TS, CXST_THETA6_TC, CXST_THETA6_A, CXST_THETA6_B);

		c_number cxst_energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);

		if(cxst_energy < (c_number) 0) {
			// derivatives called at the relevant arguments
			c_number f2D = _f2D(rstackmod, CXST_F2);
			c_number f4t1D = (use_oxDNA2_coaxial_stacking) ?
					-_f4D(t1, CXST_THETA1_T0_OXDNA2, CXST_THETA1_TS, CXST_THETA1_TC, CXST_THETA1_A, CXST_THETA1_B) - _f4D_pure_harmonic(t1, CXST_THETA1_SA, CXST_THETA1_SB) :
					-_f4D(t1, CXST_THETA1_T0_OXDNA, CXST_THETA1_TS, CXST_THETA1_TC, CXST_THETA1_A, CXST_THETA1_B) + _f4D(2 * PI - t1, CXST_THETA1_T0_OXDNA, CXST_THETA1_TS, CXST_THETA1_TC, CXST_THETA1_A, CXST_THETA1_B);
			c_number f4t4D = _f4D(t4, CXST_THETA4_T0, CXST_THETA4_TS, CXST_THETA4_TC, CXST_THETA4_A, CXST_THETA4_B);
			c_number f4t5D = _f4D(t5, CXST_THETA5_T0, CXST_THETA5_TS, CXST_THETA5_TC, CXST_THETA5_A, CXST_THETA5_B) - _f4D(PI - t5, CXST_THETA5_T0, CXST_THETA5_TS, CXST_THETA5_TC, CXST_THETA5_A, CXST_THETA5_B);
			c_number f4t6D = -_f4D(t6, CXST_THETA6_T0, CXST_THETA6_TS, CXST_THETA6_TC, CXST_THETA6_A, CXST_THETA6_B) + _f4D(PI - t6, CXST_THETA6_T0, CXST_THETA6_TS, CXST_THETA6_TC, CXST_THETA6_A, CXST_THETA6_B);

			// RADIAL PART
			Ftmp = rstackdir * (cxst_energy * f2D / f2);

			// THETA1; t1 = LRACOS (-a1 * b1);
			Ttmp -= stably_normalised(_cross(a1, b1)) * (-cxst_energy * f4t1D / f4t1);

			// TETA4; t4 = LRACOS (a3 * b3);
			Ttmp -= stably_normalised(_cross(a3, b3)) * (-cxst_energy * f4t4D / f4t4);

			// THETA5; t5 = LRACOS ( a3 * rstackdir);
			c_number part = cxst_energy * f4t5D / f4t5;
			Ftmp -= stably_normalised(a3 - rstackdir * cost5) / rstackmod * part;
			Ttmp -= stably_normalised(_cross(rstackdir, a3)) * part;

			// THETA6; t6 = LRACOS (-b3 * rstackdir);
			Ftmp -= stably_normalised(b3 + rstackdir * cost6) * (cxst_energy * f4t6D / (f4t6 * rstackmod));

			if(!use_oxDNA2_coaxial_stacking) {
				// COSPHI3
				c_number f5cosphi3D = _f5D(cosphi3, CXST_F5_PHI3);
				c_number rbackrefmodcub = CUB(rbackrefmod);

				//c_number a1b1 = a1 * b1;
				c_number a2b1 = CUDA_DOT(a2, b1);
				c_number a3b1 = CUDA_DOT(a3, b1);
				c_number ra1 = CUDA_DOT(rstackdir, a1);
				c_number ra2 = CUDA_DOT(rstackdir, a2);
				c_number ra3 = CUDA_DOT(rstackdir, a3);
				c_number rb1 = CUDA_DOT(rstackdir, b1);

				c_number parentesi = (ra3 * a2b1 - ra2 * a3b1);
				c_number dcdr = -GAMMA * parentesi * (GAMMA * (ra1 - rb1) + rstackmod) / rbackrefmodcub;
				c_number dcda1b1 = GAMMA * SQR(GAMMA) * parentesi / rbackrefmodcub;
				c_number dcda2b1 = GAMMA * ra3 / rbackrefmod;
				c_number dcda3b1 = -GAMMA * ra2 / rbackrefmod;
				c_number dcdra1 = -SQR(GAMMA) * parentesi * rstackmod / rbackrefmodcub;
				c_number dcdra2 = -GAMMA * a3b1 / rbackrefmod;
				c_number dcdra3 = GAMMA * a2b1 / rbackrefmod;
				c_number dcdrb1 = -dcdra1;

				part = cxst_energy * 2 * f5cosphi3D / f5cosphi3;

				Ftmp -= part * (rstackdir * dcdr + ((a1 - rstackdir * ra1) * dcdra1 + (a2 - rstackdir * ra2) * dcdra2 + (a3 - rstackdir * ra3) * dcdra3 + (b1 - rstackdir * rb1) * dcdrb1) / rstackmod);

				Ttmp += part * (_cross(rstackdir, a1) * dcdra1 + _cross(rstackdir, a2) * dcdra2 + _cross(rstackdir, a3) * dcdra3);

				Ttmp -= part * (_cross(a1, b1) * dcda1b1 + _cross(a2, b1) * dcda2b1 + _cross(a3, b1) * dcda3b1);
			}

			Ttmp += _cross(ppos_stack, Ftmp);
			Ftmp.w = cxst_energy;
			F += Ftmp;
		}
	}

	T += Ttmp;

	// this component stores the energy due to hydrogen bonding
	T.w = tot_Tw;
}

__global__ void dna_forces_edge_nonbonded(const c_number4 __restrict__ *poss, const GPU_quat __restrict__ *orientations, c_number4 __restrict__ *forces,
		c_number4 __restrict__ *torques, const edge_bond __restrict__ *edge_list, int n_edges, const int *is_strand_end, bool grooving,
		bool use_debye_huckel, bool use_oxDNA2_coaxial_stacking, bool update_st, CUDAStressTensor *st, const CUDABox *box) {
	if(IND >= n_edges) return;

	c_number4 dF = make_c_number4(0, 0, 0, 0);
	c_number4 dT = make_c_number4(0, 0, 0, 0);

	edge_bond b = edge_list[IND];

	// get info for particle 1
	c_number4 ppos = poss[b.from];
	bool p_is_end = (use_debye_huckel) ? is_strand_end[b.from] : false;
	// particle axes according to Allen's paper
	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[b.from], a1, a2, a3);

	// get info for particle 2
	c_number4 qpos = poss[b.to];
	bool q_is_end = (use_debye_huckel) ? is_strand_end[b.to] : false;
	c_number4 b1, b2, b3;
	get_vectors_from_quat(orientations[b.to], b1, b2, b3);

	c_number4 r = box->minimum_image(ppos, qpos);
	_particle_particle_DNA_interaction(r, ppos, a1, a2, a3, qpos, b1, b2, b3, dF, dT, grooving, use_debye_huckel, use_oxDNA2_coaxial_stacking, p_is_end, q_is_end);

	int from_index = MD_N[0] * (IND % MD_n_forces[0]) + b.from;
	int to_index = MD_N[0] * (IND % MD_n_forces[0]) + b.to;

	if(CUDA_DOT(dT, dT) > (c_number) 0.f) LR_atomicAddXYZ(&(torques[from_index]), dT);
	dT = -dT; // torque acting on particle q (which also has another contribution if dF > 0)

	if(CUDA_DOT(dF, dF) > (c_number) 0.f) {
		LR_atomicAddXYZ(&(forces[from_index]), dF);

		if(update_st) {
			CUDAStressTensor p_st;
			_update_stress_tensor<false>(p_st, r, dF);
			LR_atomicAddST(&(st[b.from]), p_st);
		}

		// Allen Eq. 6 pag 3:
		dT += _cross(r, dF);

		LR_atomicAddXYZ(&(forces[to_index]), -dF);
	}

	// the torque may be different from 0 even if dF == 0
	if(CUDA_DOT(dT, dT) > (c_number) 0.f) LR_atomicAddXYZ(&(torques[to_index]), dT);
}

// bonded interactions for edge-based approach
__global__ void dna_forces_edge_bonded(const c_number4 __restrict__ *poss, const GPU_quat __restrict__ *orientations, c_number4 __restrict__ *forces,
		c_number4 __restrict__ *torques, const LR_bonds __restrict__ *bonds, bool grooving, bool use_oxDNA2_FENE, bool use_mbf, c_number mbf_xmax,
		c_number mbf_finf, bool update_st, CUDAStressTensor *st) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = torques[IND];

	c_number4 ppos = poss[IND];
	LR_bonds bs = bonds[IND];

	CUDAStressTensor p_st;

	// particle axes according to Allen's paper
	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[IND], a1, a2, a3);

	if(bs.n3 != P_INVALID) {
		c_number4 qpos = poss[bs.n3];

		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[bs.n3], b1, b2, b3);

		c_number4 r = qpos - ppos;
		c_number4 dF = make_c_number4(0, 0, 0, 0);
		_bonded_part<true>(r, ppos, a1, a2, a3, qpos, b1, b2, b3, dF, T, grooving, use_oxDNA2_FENE, use_mbf, mbf_xmax, mbf_finf);
		_update_stress_tensor<true>(p_st, r, dF);
		F += dF;
	}
	if(bs.n5 != P_INVALID) {
		c_number4 qpos = poss[bs.n5];

		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[bs.n5], b1, b2, b3);

		c_number4 r = ppos - qpos;
		c_number4 dF = make_c_number4(0, 0, 0, 0);
		_bonded_part<false>(r, qpos, b1, b2, b3, ppos, a1, a2, a3, dF, T, grooving, use_oxDNA2_FENE, use_mbf, mbf_xmax, mbf_finf);
		_update_stress_tensor<true>(p_st, -r, dF); // -r since r here is defined as r  = ppos - qpos
		F += dF;
	}

	// we can do this because dna_forces_edge_nonbonded does not change the reference frame to the torque it calculates
	T = _vectors_transpose_c_number4_product(a1, a2, a3, T);

	forces[IND] = F;
	torques[IND] = T;

	if(update_st) {
		st[IND] += p_st;
	}
}

__global__ void dna_forces(const c_number4 __restrict__ *poss, const GPU_quat __restrict__ *orientations, c_number4 __restrict__ *forces,
		c_number4 __restrict__ *torques, const int *matrix_neighs,	const int *number_neighs, const LR_bonds __restrict__ *bonds,
		bool grooving, bool use_debye_huckel, bool use_oxDNA2_coaxial_stacking, bool use_oxDNA2_FENE, bool use_mbf, c_number mbf_xmax,
		c_number mbf_finf, bool update_st, CUDAStressTensor *st, const CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = make_c_number4(0, 0, 0, 0);
	c_number4 ppos = poss[IND];
	LR_bonds pbonds = bonds[IND];
	bool p_is_end = (pbonds.n3 == P_INVALID || pbonds.n5 == P_INVALID);

	CUDAStressTensor p_st;

	// particle axes according to Allen's paper
	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[IND], a1, a2, a3);

	if(pbonds.n3 != P_INVALID) {
		c_number4 qpos = poss[pbonds.n3];
		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[pbonds.n3], b1, b2, b3);

		c_number4 r = qpos - ppos;
		c_number4 dF = make_c_number4(0, 0, 0, 0);
		_bonded_part<true>(r, ppos, a1, a2, a3, qpos, b1, b2, b3, dF, T, grooving, use_oxDNA2_FENE, use_mbf, mbf_xmax, mbf_finf);
		_update_stress_tensor<true>(p_st, r, dF);
		F += dF;
	}
	if(pbonds.n5 != P_INVALID) {
		c_number4 qpos = poss[pbonds.n5];
		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[pbonds.n5], b1, b2, b3);

		c_number4 r = ppos - qpos;
		c_number4 dF = make_c_number4(0, 0, 0, 0);
		_bonded_part<false>(r, qpos, b1, b2, b3, ppos, a1, a2, a3, dF, T, grooving, use_oxDNA2_FENE, use_mbf, mbf_xmax, mbf_finf);
		_update_stress_tensor<true>(p_st, -r, dF);
		F += dF;
	}

	const int type = get_particle_type(ppos);
	int num_neighs = NUMBER_NEIGHBOURS(IND, number_neighs);

	T.w = (c_number) 0;
	for(int j = 0; j < num_neighs; j++) {
		int k_index = NEXT_NEIGHBOUR(IND, j, matrix_neighs);

		if(k_index != IND && pbonds.n3 != k_index && pbonds.n5 != k_index) {
			c_number4 qpos = poss[k_index];
			c_number4 r = box->minimum_image(ppos, qpos);

			c_number4 b1, b2, b3;
			get_vectors_from_quat(orientations[k_index], b1, b2, b3);
			LR_bonds qbonds = bonds[k_index];
			bool q_is_end = (qbonds.n3 == P_INVALID || qbonds.n5 == P_INVALID);

			c_number4 dF = make_c_number4(0, 0, 0, 0);
			_particle_particle_DNA_interaction(r, ppos, a1, a2, a3, qpos, b1, b2, b3, dF, T, grooving, use_debye_huckel, use_oxDNA2_coaxial_stacking, p_is_end, q_is_end);
			_update_stress_tensor<true>(p_st, r, dF);
			F += dF;
		}
	}

	T = _vectors_transpose_c_number4_product(a1, a2, a3, T);

	forces[IND] = F;
	torques[IND] = T;

	if(update_st) {
		st[IND] = p_st;
	}
}

//FFS order parameter pre-calculations

// check whether a particular pair of particles have hydrogen bonding energy lower than a given threshold hb_threshold (which may vary)

__global__ void hb_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, float *hb_energies, int n_threads, bool *region_is_nearhb, CUDABox *box) {
	if(IND >= n_threads) return;

	int pind = op_pairs1[IND];
	int qind = op_pairs2[IND];
	// get distance between this nucleotide pair's "com"s
	c_number4 ppos = poss[pind];
	c_number4 qpos = poss[qind];
	c_number4 r = box->minimum_image(ppos, qpos);

	// check whether hb energy is below a certain threshold for this nucleotide pair
	int ptype = get_particle_type(ppos);
	int qtype = get_particle_type(qpos);
	int pbtype = get_particle_btype(ppos);
	int qbtype = get_particle_btype(qpos);
	int int_type = pbtype + qbtype;

	GPU_quat po = orientations[pind];
	GPU_quat qo = orientations[qind];

	//This gets an extra two vectors that are not needed, but the function doesn't seem to be called at all, so should make no difference. 
	c_number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);

	get_vectors_from_quat(qo, b1, b2, b3);

	c_number4 ppos_base = POS_BASE * a1;
	c_number4 qpos_base = POS_BASE * b1;

	// HYDROGEN BONDING
	c_number hb_energy = (c_number) 0;
	c_number4 rhydro = r + qpos_base - ppos_base;
	c_number rhydromodsqr = CUDA_DOT(rhydro, rhydro);
	if(int_type == 3 && SQR(HYDR_RCLOW) < rhydromodsqr && rhydromodsqr < SQR(HYDR_RCHIGH)) {
		// versor and magnitude of the base-base separation
		c_number rhydromod = sqrtf(rhydromodsqr);
		c_number4 rhydrodir = rhydro / rhydromod;

		// angles involved in the HB interaction
		c_number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));
		c_number t2 = CUDA_LRACOS(-CUDA_DOT(b1, rhydrodir));
		c_number t3 = CUDA_LRACOS(CUDA_DOT(a1, rhydrodir));
		c_number t4 = CUDA_LRACOS(CUDA_DOT(a3, b3));
		c_number t7 = CUDA_LRACOS(-CUDA_DOT(rhydrodir, b3));
		c_number t8 = CUDA_LRACOS(CUDA_DOT(rhydrodir, a3));

		// functions called at their relevant arguments
		c_number f1 = _f1(rhydromod, HYDR_F1, ptype, qtype);
		c_number f4t1 = _f4(t1, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC, HYDR_THETA1_A, HYDR_THETA1_B);
		c_number f4t2 = _f4(t2, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC, HYDR_THETA2_A, HYDR_THETA2_B);
		c_number f4t3 = _f4(t3, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC, HYDR_THETA3_A, HYDR_THETA3_B);
		c_number f4t4 = _f4(t4, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC, HYDR_THETA4_A, HYDR_THETA4_B);
		c_number f4t7 = _f4(t7, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC, HYDR_THETA7_A, HYDR_THETA7_B);
		c_number f4t8 = _f4(t8, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC, HYDR_THETA8_A, HYDR_THETA8_B);

		hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
	}
	// END HYDROGEN BONDING

	hb_energies[IND] = hb_energy;
}

__global__ void near_hb_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, bool *nearly_bonded_array, int n_threads, bool *region_is_nearhb, CUDABox *box) {
	if(IND >= n_threads) return;

	int pind = op_pairs1[IND];
	int qind = op_pairs2[IND];
	// get distance between this nucleotide pair's "com"s
	c_number4 ppos = poss[pind];
	c_number4 qpos = poss[qind];
	c_number4 r = box->minimum_image(ppos, qpos);

	// check whether hb energy is below a certain threshold for this nucleotide pair
	int ptype = get_particle_type(ppos);
	int qtype = get_particle_type(qpos);
	int pbtype = get_particle_btype(ppos);
	int qbtype = get_particle_btype(qpos);
	int int_type = pbtype + qbtype;

	GPU_quat po = orientations[pind];
	GPU_quat qo = orientations[qind];

	//This gets extra a2 and b2 vectors that aren't needed. get_vectors_from_quat could easily be modified to only return the relevant vectors, but it will make computationally very little difference, since most of the same c_numbers need to be calculated anyway. Perhaps worth changing for memory considerations, however
	c_number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);
	get_vectors_from_quat(qo, b1, b2, b3);

	c_number4 ppos_base = POS_BASE * a1;
	c_number4 qpos_base = POS_BASE * b1;

	// HYDROGEN BONDING
	c_number4 rhydro = r + qpos_base - ppos_base;
	c_number rhydromodsqr = CUDA_DOT(rhydro, rhydro);

	int total_nonzero;
	bool nearly_bonded = false;
	if(int_type == 3 && SQR(HYDR_RCLOW) < rhydromodsqr && rhydromodsqr < SQR(HYDR_RCHIGH)) {
		// versor and magnitude of the base-base separation
		c_number rhydromod = sqrtf(rhydromodsqr);
		c_number4 rhydrodir = rhydro / rhydromod;

		// angles involved in the HB interaction
		c_number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));
		c_number t2 = CUDA_LRACOS(-CUDA_DOT(b1, rhydrodir));
		c_number t3 = CUDA_LRACOS(CUDA_DOT(a1, rhydrodir));
		c_number t4 = CUDA_LRACOS(CUDA_DOT(a3, b3));
		c_number t7 = CUDA_LRACOS(-CUDA_DOT(rhydrodir, b3));
		c_number t8 = CUDA_LRACOS(CUDA_DOT(rhydrodir, a3));

		// functions called at their relevant arguments
		c_number f1 = _f1(rhydromod, HYDR_F1, ptype, qtype);
		c_number f4t1 = _f4(t1, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC, HYDR_THETA1_A, HYDR_THETA1_B);
		c_number f4t2 = _f4(t2, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC, HYDR_THETA2_A, HYDR_THETA2_B);
		c_number f4t3 = _f4(t3, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC, HYDR_THETA3_A, HYDR_THETA3_B);
		c_number f4t4 = _f4(t4, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC, HYDR_THETA4_A, HYDR_THETA4_B);
		c_number f4t7 = _f4(t7, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC, HYDR_THETA7_A, HYDR_THETA7_B);
		c_number f4t8 = _f4(t8, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC, HYDR_THETA8_A, HYDR_THETA8_B);

		// the nearly bonded order parameter requires either all or all but one of these factors to be non-zero
		bool f1not0 = f1 < 0 ? 1 : 0;
		bool f4t1not0 = f4t1 > 0 ? 1 : 0;
		bool f4t2not0 = f4t2 > 0 ? 1 : 0;
		bool f4t3not0 = f4t3 > 0 ? 1 : 0;
		bool f4t4not0 = f4t4 > 0 ? 1 : 0;
		bool f4t7not0 = f4t7 > 0 ? 1 : 0;
		bool f4t8not0 = f4t8 > 0 ? 1 : 0;
		total_nonzero = f1not0 + f4t1not0 + f4t2not0 + f4t3not0 + f4t4not0 + f4t7not0 + f4t8not0;
		if(total_nonzero >= 6) nearly_bonded = true;
	}
	// END HYDROGEN BONDING

	nearly_bonded_array[IND] = nearly_bonded;
}

// compute the distance between a pair of particles
__global__ void dist_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, c_number *op_dists, int n_threads, CUDABox *box) {
	if(IND >= n_threads) return;

	int pind = op_pairs1[IND];
	int qind = op_pairs2[IND];

	// get distance between this nucleotide pair's "com"s
	c_number4 ppos = poss[pind];
	c_number4 qpos = poss[qind];
	c_number4 r = box->minimum_image(ppos, qpos);

	GPU_quat po = orientations[pind];
	GPU_quat qo = orientations[qind];

	//This gets extra a2 and b2 vectors that aren't needed. get_vectors_from_quat could easily be modified to only return the relevant vectors, but it will make computationally very little difference, since most of the same c_numbers need to be calculated anyway. Perhaps worth changing for memory considerations, however 
	c_number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);
	get_vectors_from_quat(qo, b1, b2, b3);

	c_number4 ppos_base = POS_BASE * a1;
	c_number4 qpos_base = POS_BASE * b1;

	c_number4 rbase = r + qpos_base - ppos_base;
	op_dists[IND] = _module(rbase);
}

__global__ void init_DNA_strand_ends(int *is_strand_end, const LR_bonds __restrict__ *bonds, int N) {
	if(IND >= N) return;

	LR_bonds pbonds = bonds[IND];
	is_strand_end[IND] = (pbonds.n3 == P_INVALID || pbonds.n5 == P_INVALID);
}

