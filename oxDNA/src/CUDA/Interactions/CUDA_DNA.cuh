/* System constants */
__constant__ int MD_N[1];
__constant__ float MD_box_side[1];
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

template<typename number, typename number4>
__forceinline__ __device__ void _excluded_volume(const number4 &r, number4 &F, number sigma, number rstar, number b, number rc) {
	number rsqr = CUDA_DOT(r, r);

	F.x = F.y = F.z = F.w = (number) 0.f;
	if(rsqr < SQR(rc)) {
		if(rsqr > SQR(rstar)) {
			number rmod = sqrt(rsqr);
			number rrc = rmod - rc;
			number fmod = 2.f * EXCL_EPS * b * rrc / rmod;
			F.x = r.x * fmod;
			F.y = r.y * fmod;
			F.z = r.z * fmod;
			F.w = EXCL_EPS * b * SQR(rrc);
		}
		else {
			number lj_part = CUB(SQR(sigma)/rsqr);
			number fmod = 24.f * EXCL_EPS * (lj_part - 2.f*SQR(lj_part)) / rsqr;
			F.x = r.x * fmod;
			F.y = r.y * fmod;
			F.z = r.z * fmod;
			F.w = 4.f * EXCL_EPS * (SQR(lj_part) - lj_part);
		}
	}
}

template<typename number>
__forceinline__ __device__ number _f1(number r, int type, int n3, int n5) {
	number val = (number) 0.f;
	if(r < MD_F1_RCHIGH[type]) {
		int eps_index = 25 * type + n3 * 5 + n5;
		if(r > MD_F1_RHIGH[type]) {
			val = MD_F1_EPS[eps_index] * MD_F1_BHIGH[type] * SQR(r - MD_F1_RCHIGH[type]);
		}
		else if(r > MD_F1_RLOW[type]) {
			number tmp = 1.f - expf(-(r - MD_F1_R0[type]) * MD_F1_A[type]);
			val = MD_F1_EPS[eps_index] * SQR(tmp) - MD_F1_SHIFT[eps_index];
		}
		else if(r > MD_F1_RCLOW[type]) {
			val = MD_F1_EPS[eps_index] * MD_F1_BLOW[type] * SQR(r - MD_F1_RCLOW[type]);
		}
	}

	return val;
}

template<typename number>
__forceinline__ __device__ number _f1D(number r, int type, int n3, int n5) {
	number val = (number) 0.f;
	int eps_index = 0;
	if(r < MD_F1_RCHIGH[type]) {
		eps_index = 25 * type + n3 * 5 + n5;
		if(r > MD_F1_RHIGH[type]) {
			val = 2.f * MD_F1_BHIGH[type] * (r - MD_F1_RCHIGH[type]);
		}
		else if(r > MD_F1_RLOW[type]) {
			number tmp = expf(-(r - MD_F1_R0[type]) * MD_F1_A[type]);
			val = 2.f * (1.f - tmp) * tmp * MD_F1_A[type];
		}
		else if(r > MD_F1_RCLOW[type]) {
			val = 2.f * MD_F1_BLOW[type] * (r - MD_F1_RCLOW[type]);
		}
	}

	return MD_F1_EPS[eps_index] * val;
}

template<typename number>
__forceinline__ __device__ number _f2(number r, int type) {
    number val = (number) 0.f;
    if (r < MD_F2_RCHIGH[type]) {
	    if (r > MD_F2_RHIGH[type]) {
		    val = MD_F2_K[type] * MD_F2_BHIGH[type] * SQR(r - MD_F2_RCHIGH[type]);
	    }
	    else if (r > MD_F2_RLOW[type]) {
		    val = (MD_F2_K[type] * 0.5f) * (SQR(r - MD_F2_R0[type]) - SQR(MD_F2_RC[type] - MD_F2_R0[type]));
	    }
	    else if (r > MD_F2_RCLOW[type]) {
		    val = MD_F2_K[type] * MD_F2_BLOW[type] * SQR(r - MD_F2_RCLOW[type]);
	    }
    }
    return val;
}

template<typename number>
__forceinline__ __device__ number _f2D(number r, int type) {
    number val = (number) 0.f;
    if (r < MD_F2_RCHIGH[type]) {
	    if (r > MD_F2_RHIGH[type]) {
		    val = 2.f * MD_F2_K[type] * MD_F2_BHIGH[type] * (r - MD_F2_RCHIGH[type]);
	    }
	    else if (r > MD_F2_RLOW[type]) {
		    val = MD_F2_K[type] * (r - MD_F2_R0[type]);
	    }
	    else if (r > MD_F2_RCLOW[type]) {
		    val = 2.f * MD_F2_K[type] * MD_F2_BLOW[type] * (r - MD_F2_RCLOW[type]);
	    }
    }
    return val;
}

template<typename number>
__forceinline__ __device__ number _f4(number t, float t0, float ts, float tc, float a, float b) {
	number val = (number) 0.f;
	t -= t0;
	if(t < 0) t = -t;

	if(t < tc) {
		if(t > ts) {
			// smoothing
			val = b * SQR(tc - t);
		}
		else val = (number) 1.f - a * SQR(t);
	}

	return val;
}

template<typename number>
__forceinline__ __device__ number _f4_pure_harmonic(number t, float a, float b) {
	// for getting a f4t1 function with a continuous derivative that is less disruptive to the potential
	number val = (number) 0.f;
	t -= b;
	if(t < 0) val = (number) 0.f;
	else val = (number) a * SQR(t);

	return val;
}

template<typename number>
__forceinline__ __device__ number _f4Dsin(number t, float t0, float ts, float tc, float a, float b) {
	number val = (number) 0.f;
	number tt0 = t - t0;
	// this function is a parabola centered in t0. If tt0 < 0 then the value of the function
	// is the same but the value of its derivative has the opposite sign, so m = -1
	number m = copysignf((number)1.f, tt0);
	tt0 = copysignf(tt0, (number)1.f);

	if(tt0 < tc) {
		number sint = sinf(t);
		if(tt0 > ts) {
			// smoothing
			val = b * (tt0 - tc) / sint;
		}
		else {
			if(SQR(sint) > 1e-12f) val = -a * tt0 / sint;
			else val = -a;
		}
	}

	return 2.f * m * val;
}

template<typename number>
__forceinline__ __device__ number _f4Dsin_pure_harmonic(number t, float a, float b) {
	// for getting a f4t1 function with a continuous derivative that is less disruptive to the potential
	number val = (number) 0.f;
	number tt0 = t - b;
	if(tt0 < 0) val = (number) 0.f;
	else {
		number sint = sin(t);
		if (SQR(sint) > 1e-12) val = (number) 2 * a * tt0 / sint;
		else val = (number) 2 * a;
	}
	
	return val;
}

template<typename number>
__forceinline__ __device__ number _f5(number f, int type) {
	number val = (number) 0.f;

	if(f > MD_F5_PHI_XC[type]) {
		if(f < MD_F5_PHI_XS[type]) {
			val = MD_F5_PHI_B[type] * SQR(MD_F5_PHI_XC[type] - f);
		}
		else if(f < 0.f) {
			val = (number) 1.f - MD_F5_PHI_A[type] * SQR(f);
		}
		else val = 1.f;
	}

	return val;
}

template<typename number>
__forceinline__ __device__ number _f5D(number f, int type) {
	number val = (number) 0.f;

	if(f > MD_F5_PHI_XC[type]) {
		if(f < MD_F5_PHI_XS[type]) {
			val = 2.f * MD_F5_PHI_B[type] * (f - MD_F5_PHI_XC[type]);
		}
		else if(f < 0.f) {
			val = (number) -2.f * MD_F5_PHI_A[type] * f;
		}
	}

	return val;
}

template <typename number, typename number4>
__device__ number4 minimum_image(const number4 &r_i, const number4 &r_j) {
	number dx = r_j.x - r_i.x;
	number dy = r_j.y - r_i.y;
	number dz = r_j.z - r_i.z;

	dx -= floorf(dx/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];
	dy -= floorf(dy/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];
	dz -= floorf(dz/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];

	return make_number4<number, number4>(dx, dy, dz, (number) 0.f);
}

template <typename number, typename number4>
__device__ number quad_minimum_image_dist(const number4 &r_i, const number4 &r_j) {
	number dx = r_j.x - r_i.x;
	number dy = r_j.y - r_i.y;
	number dz = r_j.z - r_i.z;

	dx -= floorf(dx/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];
	dy -= floorf(dy/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];
	dz -= floorf(dz/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];

	return dx*dx + dy*dy + dz*dz;
}

template <typename number, typename number4, bool qIsN3>
__device__ void _bonded_excluded_volume(number4 &r, number4 &n3pos_base, number4 &n3pos_back, number4 &n5pos_base, number4 &n5pos_back, number4 &F, number4 &T) {
	number4 Ftmp;
	// BASE-BASE
	number4 rcenter = r + n3pos_base - n5pos_base;
	_excluded_volume(rcenter, Ftmp, EXCL_S2, EXCL_R2, EXCL_B2, EXCL_RC2);
	number4 torquep1 = (qIsN3) ? _cross<number, number4>(n5pos_base, Ftmp) : _cross<number, number4>(n3pos_base, Ftmp);
	F += Ftmp;

	// n5-BASE vs. n3-BACK
	rcenter = r + n3pos_back - n5pos_base;
	_excluded_volume(rcenter, Ftmp, EXCL_S3, EXCL_R3, EXCL_B3, EXCL_RC3);
	number4 torquep2 = (qIsN3) ? _cross<number, number4>(n5pos_base, Ftmp) : _cross<number, number4>(n3pos_back, Ftmp);
	F += Ftmp;

	// n5-BACK vs. n3-BASE
	rcenter = r + n3pos_base - n5pos_back;
	_excluded_volume(rcenter, Ftmp, EXCL_S4, EXCL_R4, EXCL_B4, EXCL_RC4);
	number4 torquep3 = (qIsN3) ? _cross<number, number4>(n5pos_back, Ftmp) : _cross<number, number4>(n3pos_base, Ftmp);
	F += Ftmp;

	T += torquep1 + torquep2 + torquep3;
}

template <typename number, typename number4, bool qIsN3>
__device__ void _bonded_part(number4 &n5pos, number4 &n5x, number4 &n5y, number4 &n5z,number4 &n3pos,
			     number4 &n3x, number4 &n3y, number4 &n3z, number4 &F, number4 &T, bool grooving, bool use_oxDNA2_FENE) {

	int n3type = get_particle_type<number, number4>(n3pos);
	int n5type = get_particle_type<number, number4>(n5pos);

	number4 r = make_number4<number, number4>(n3pos.x - n5pos.x, n3pos.y - n5pos.y, n3pos.z - n5pos.z, (number) 0);

	number4 n5pos_back;
	if(grooving) n5pos_back = n5x * POS_MM_BACK1 + n5y * POS_MM_BACK2;
	else n5pos_back = n5x * POS_BACK;

	number4 n5pos_base = n5x * POS_BASE;
	number4 n5pos_stack = n5x * POS_STACK;

	number4 n3pos_back;
	if(grooving) n3pos_back = n3x * POS_MM_BACK1 + n3y * POS_MM_BACK2;
	else n3pos_back = n3x * POS_BACK;
	number4 n3pos_base = n3x * POS_BASE;
	number4 n3pos_stack = n3x * POS_STACK;

	number4 rback = r + n3pos_back - n5pos_back;
	number rbackmod = _module<number, number4>(rback);
	number rbackr0;
	if (use_oxDNA2_FENE) rbackr0 = rbackmod - FENE_R0_OXDNA2;
	else rbackr0 = rbackmod - FENE_R0_OXDNA;

	number4 Ftmp = rback * ((FENE_EPS * rbackr0  / (FENE_DELTA2 - SQR(rbackr0))) / rbackmod);
	Ftmp.w = -FENE_EPS * ((number)0.5f) * logf(1 - SQR(rbackr0) / FENE_DELTA2);

	number4 Ttmp = (qIsN3) ? _cross<number, number4>(n5pos_back, Ftmp) : _cross<number, number4>(n3pos_back, Ftmp);
	// EXCLUDED VOLUME
	_bonded_excluded_volume<number, number4, qIsN3>(r, n3pos_base, n3pos_back, n5pos_base, n5pos_back, Ftmp, Ttmp);

	if(qIsN3) {
		F += Ftmp;
		T += Ttmp;
	}
	else {
		F -= Ftmp;
		T -= Ttmp;
	}

	// STACKING
	number4 rstack = r + n3pos_stack - n5pos_stack;
	number rstackmod = _module<number, number4>(rstack);
	number4 rstackdir = make_number4<number, number4>(rstack.x / rstackmod, rstack.y / rstackmod, rstack.z / rstackmod, 0);
	// This is the position the backbone would have with major-minor grooves the same width.
	// We need to do this to implement different major-minor groove widths because rback is
	// used as a reference point for things that have nothing to do with the actual backbone
	// position (in this case, the stacking interaction).
	number4 rbackref = r + n3x * POS_BACK - n5x * POS_BACK;
	number rbackrefmod = _module<number, number4>(rbackref);

	number t4 = CUDA_LRACOS(CUDA_DOT(n3z, n5z));
	number t5 = CUDA_LRACOS(CUDA_DOT(n5z, rstackdir));
	number t6 = CUDA_LRACOS(-CUDA_DOT(n3z, rstackdir));
	number cosphi1 = CUDA_DOT(n5y, rbackref) / rbackrefmod;
	number cosphi2 = CUDA_DOT(n3y, rbackref) / rbackrefmod;

	// functions
	number f1 = _f1(rstackmod, STCK_F1, n3type, n5type);
	number f4t4 = _f4(t4, STCK_THETA4_T0, STCK_THETA4_TS, STCK_THETA4_TC, STCK_THETA4_A, STCK_THETA4_B);
	number f4t5 = _f4(PI - t5, STCK_THETA5_T0, STCK_THETA5_TS, STCK_THETA5_TC, STCK_THETA5_A, STCK_THETA5_B);
	number f4t6 = _f4(t6, STCK_THETA6_T0, STCK_THETA6_TS, STCK_THETA6_TC, STCK_THETA6_A, STCK_THETA6_B);
	number f5phi1 = _f5(cosphi1, STCK_F5_PHI1);
	number f5phi2 = _f5(cosphi2, STCK_F5_PHI2);

	number energy = f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;

	if(energy != (number) 0) {
		// and their derivatives
		number f1D = _f1D(rstackmod, STCK_F1, n3type, n5type);
		number f4t4Dsin = _f4Dsin(t4, STCK_THETA4_T0, STCK_THETA4_TS, STCK_THETA4_TC, STCK_THETA4_A, STCK_THETA4_B);
			number f4t5Dsin = _f4Dsin(PI - t5, STCK_THETA5_T0, STCK_THETA5_TS, STCK_THETA5_TC, STCK_THETA5_A, STCK_THETA5_B);
			number f4t6Dsin = _f4Dsin(t6, STCK_THETA6_T0, STCK_THETA6_TS, STCK_THETA6_TC, STCK_THETA6_A, STCK_THETA6_B);
		number f5phi1D = _f5D(cosphi1, STCK_F5_PHI1);
		number f5phi2D = _f5D(cosphi2, STCK_F5_PHI2);

		// RADIAL
		Ftmp = rstackdir * (energy * f1D / f1);

		// THETA 5
		Ftmp += (n5z - cosf(t5) * rstackdir) * (energy * f4t5Dsin / (f4t5 * rstackmod));

		// THETA 6
		Ftmp += (n3z + cosf(t6) * rstackdir) * (energy * f4t6Dsin / (f4t6 * rstackmod));

		// COS PHI 1
		// here particle p is referred to using the a while particle q is referred with the b
		number ra2 = CUDA_DOT(rstackdir, n5y);
		number ra1 = CUDA_DOT(rstackdir, n5x);
		number rb1 = CUDA_DOT(rstackdir, n3x);
		number a2b1 = CUDA_DOT(n5y, n3x);
		number dcosphi1dr = (SQR(rstackmod)*ra2 - ra2*SQR(rbackrefmod) - rstackmod*(a2b1 + ra2*(-ra1 + rb1))*GAMMA + a2b1*(-ra1 + rb1)*SQR(GAMMA))/(SQR(rbackrefmod)*rbackrefmod);
		number dcosphi1dra1 = rstackmod*GAMMA*(rstackmod*ra2 - a2b1*GAMMA)/(SQR(rbackrefmod)*rbackrefmod);
		number dcosphi1dra2 = -rstackmod / rbackrefmod;
		number dcosphi1drb1 = -(rstackmod*GAMMA*(rstackmod*ra2 - a2b1*GAMMA))/(SQR(rbackrefmod)*rbackrefmod);
		number dcosphi1da1b1 = SQR(GAMMA)*(-rstackmod*ra2 + a2b1*GAMMA)/(SQR(rbackrefmod)*rbackrefmod);
		number dcosphi1da2b1 = GAMMA / rbackrefmod;

		number force_part_phi1 = energy * f5phi1D / f5phi1;

		Ftmp -= (rstackdir * dcosphi1dr +
			    ((n5y - ra2*rstackdir) * dcosphi1dra2 +
				(n5x - ra1*rstackdir) * dcosphi1dra1 +
				(n3x - rb1*rstackdir) * dcosphi1drb1) / rstackmod) * force_part_phi1;

		// COS PHI 2
		// here particle p -> b, particle q -> a
		ra2 = CUDA_DOT(rstackdir, n3y);
		ra1 = rb1;
		rb1 = CUDA_DOT(rstackdir, n5x);
		a2b1 = CUDA_DOT(n3y, n5x);
		number dcosphi2dr = ((rstackmod*ra2 + a2b1*GAMMA)*(rstackmod + (rb1 - ra1)*GAMMA) - ra2*SQR(rbackrefmod))/(SQR(rbackrefmod)*rbackrefmod);
		number dcosphi2dra1 = -rstackmod*GAMMA*(rstackmod*ra2 + a2b1*GAMMA)/(SQR(rbackrefmod)*rbackrefmod);
		number dcosphi2dra2 = -rstackmod / rbackrefmod;
		number dcosphi2drb1 = (rstackmod*GAMMA*(rstackmod*ra2 + a2b1*GAMMA))/(SQR(rbackrefmod)*rbackrefmod);
		number dcosphi2da1b1 = -SQR(GAMMA)*(rstackmod*ra2 + a2b1*GAMMA)/(SQR(rbackrefmod)*rbackrefmod);
		number dcosphi2da2b1 = -GAMMA / rbackrefmod;

		number force_part_phi2 = energy * f5phi2D / f5phi2;

		Ftmp -= (rstackdir * dcosphi2dr +
			    ((n3y - rstackdir * ra2) * dcosphi2dra2 +
				(n3x - rstackdir * ra1) * dcosphi2dra1 +
				(n5x - rstackdir * rb1) * dcosphi2drb1) / rstackmod) * force_part_phi2;

		if(qIsN3) Ttmp = _cross<number, number4>(n5pos_stack, Ftmp);
		else Ttmp = _cross<number, number4>(n3pos_stack, Ftmp);

		// THETA 4
		Ttmp += _cross<number, number4>(n3z, n5z) * (-energy * f4t4Dsin / f4t4);

		// PHI 1 & PHI 2
		if(qIsN3) {
			Ttmp += (-force_part_phi1 * dcosphi1dra2) * _cross<number, number4>(rstackdir, n5y)
				-_cross<number, number4>(rstackdir, n5x) * force_part_phi1 * dcosphi1dra1;

			Ttmp += (-force_part_phi2 * dcosphi2drb1) * _cross<number, number4>(rstackdir, n5x);
		}
		else {
			Ttmp += force_part_phi1 * dcosphi1drb1 * _cross<number, number4>(rstackdir, n3x);

			Ttmp += force_part_phi2 * dcosphi2dra2 * _cross<number, number4>(rstackdir, n3y) +
				force_part_phi2 * dcosphi2dra1 * _cross<number, number4>(rstackdir, n3x);
		}

		Ttmp += force_part_phi1 * dcosphi1da2b1 * _cross<number, number4>(n5y, n3x)
			+ _cross<number, number4>(n5x, n3x) * force_part_phi1 * dcosphi1da1b1;

		Ttmp += force_part_phi2 * dcosphi2da2b1 * _cross<number, number4>(n5x, n3y) +
			_cross<number, number4>(n5x, n3x) * force_part_phi2 * dcosphi2da1b1;

		Ftmp.w = energy;
		if(qIsN3) {
			// THETA 5
			Ttmp += _cross<number, number4>(rstackdir, n5z) * energy * f4t5Dsin / f4t5;

			T += Ttmp;
			F += Ftmp;
		}
		else {
			// THETA 6
			Ttmp += _cross<number, number4>(rstackdir, n3z) * (-energy * f4t6Dsin / f4t6);

			T -= Ttmp;
			F -= Ftmp;
		}
	}
}

template <typename number, typename number4>
__device__ void _particle_particle_interaction(number4 ppos, number4 a1, number4 a2, number4 a3, number4 qpos, number4 b1, number4 b2, number4 b3, number4 &F, number4 &T, bool grooving, bool use_debye_huckel, bool use_oxDNA2_coaxial_stacking, LR_bonds pbonds, LR_bonds qbonds, int pind, int qind) {
	int ptype = get_particle_type<number, number4>(ppos);
	int qtype = get_particle_type<number, number4>(qpos);
	int pbtype = get_particle_btype<number, number4>(ppos);
	int qbtype = get_particle_btype<number, number4>(qpos);
	int int_type = pbtype + qbtype;

	number4 r = minimum_image<number, number4>(ppos, qpos);

	number4 ppos_back;
	if(grooving) ppos_back = POS_MM_BACK1 * a1 + POS_MM_BACK2 * a2;
	else ppos_back = POS_BACK * a1;
	number4 ppos_base = POS_BASE * a1;
	number4 ppos_stack = POS_STACK * a1;

	number4 qpos_back;
	if(grooving) qpos_back = POS_MM_BACK1 * b1 + POS_MM_BACK2 * b2;
	else qpos_back = POS_BACK * b1;
	number4 qpos_base = POS_BASE * b1;
	number4 qpos_stack = POS_STACK * b1;

	number old_Tw = T.w;

	// excluded volume
	// BACK-BACK
	number4 Ftmp = make_number4<number, number4>(0, 0, 0, 0);
	number4 rbackbone = r + qpos_back - ppos_back;
	_excluded_volume(rbackbone, Ftmp, EXCL_S1, EXCL_R1, EXCL_B1, EXCL_RC1);
	number4 Ttmp = _cross<number, number4>(ppos_back, Ftmp);
	_bonded_excluded_volume<number, number4, true>(r, qpos_base, qpos_back, ppos_base, ppos_back, Ftmp, Ttmp);

	F += Ftmp;

	// HYDROGEN BONDING
	number hb_energy = (number) 0;
	number4 rhydro = r + qpos_base - ppos_base;
	number rhydromodsqr = CUDA_DOT(rhydro, rhydro);
	if(int_type == 3 && SQR(HYDR_RCLOW) < rhydromodsqr && rhydromodsqr < SQR(HYDR_RCHIGH)) {
		number hb_multi = (abs(qbtype) >= 300 && abs(pbtype) >= 300) ? MD_hb_multi[0] : 1.f;
		// versor and magnitude of the base-base separation
	  	number rhydromod = sqrtf(rhydromodsqr);
	  	number4 rhydrodir = rhydro / rhydromod;

		// angles involved in the HB interaction
		number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));
		number t2 = CUDA_LRACOS(-CUDA_DOT(b1, rhydrodir));
		number t3 = CUDA_LRACOS(CUDA_DOT(a1, rhydrodir));
		number t4 = CUDA_LRACOS(CUDA_DOT(a3, b3));
		number t7 = CUDA_LRACOS(-CUDA_DOT(rhydrodir, b3));
		number t8 = CUDA_LRACOS(CUDA_DOT(rhydrodir, a3));

	 	 // functions called at their relevant arguments
		number f1 = hb_multi * _f1(rhydromod, HYDR_F1, ptype, qtype);
		number f4t1 = _f4(t1, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC, HYDR_THETA1_A, HYDR_THETA1_B);
		number f4t2 = _f4(t2, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC, HYDR_THETA2_A, HYDR_THETA2_B);
		number f4t3 = _f4(t3, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC, HYDR_THETA3_A, HYDR_THETA3_B);
		number f4t4 = _f4(t4, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC, HYDR_THETA4_A, HYDR_THETA4_B);
		number f4t7 = _f4(t7, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC, HYDR_THETA7_A, HYDR_THETA7_B);
		number f4t8 = _f4(t8, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC, HYDR_THETA8_A, HYDR_THETA8_B);

		hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

		if(hb_energy < (number) 0) {
			// derivatives called at the relevant arguments
			number f1D = hb_multi * _f1D(rhydromod, HYDR_F1, ptype, qtype);
			number f4t1Dsin = - _f4Dsin(t1, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC, HYDR_THETA1_A, HYDR_THETA1_B);
			number f4t2Dsin = - _f4Dsin(t2, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC, HYDR_THETA2_A, HYDR_THETA2_B);
			number f4t3Dsin = _f4Dsin(t3, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC, HYDR_THETA3_A, HYDR_THETA3_B);
			number f4t4Dsin = _f4Dsin(t4, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC, HYDR_THETA4_A, HYDR_THETA4_B);
			number f4t7Dsin = - _f4Dsin(t7, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC, HYDR_THETA7_A, HYDR_THETA7_B);
			number f4t8Dsin = _f4Dsin(t8, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC, HYDR_THETA8_A, HYDR_THETA8_B);

			// RADIAL PART
			Ftmp = rhydrodir * hb_energy * f1D / f1;

			// TETA4; t4 = LRACOS (a3 * b3);
			Ttmp -= _cross<number, number4>(a3, b3) * (-hb_energy * f4t4Dsin / f4t4);

			// TETA1; t1 = LRACOS (-a1 * b1);
			Ttmp -= _cross<number, number4>(a1, b1) * (- hb_energy * f4t1Dsin / f4t1);

			// TETA2; t2 = LRACOS (-b1 * rhydrodir);
			Ftmp -= (b1 + rhydrodir * cosf(t2)) * (hb_energy * f4t2Dsin / (f4t2 * rhydromod));

			// TETA3; t3 = LRACOS (a1 * rhydrodir);
			number part = - hb_energy * f4t3Dsin / f4t3;
			Ftmp -= (a1 - rhydrodir * cosf(t3)) * (-part / rhydromod);
			Ttmp += _cross<number, number4>(rhydrodir, a1) * part;

			// THETA7; t7 = LRACOS (-rhydrodir * b3);
			Ftmp -= (b3 + rhydrodir * cosf(t7)) * (hb_energy * f4t7Dsin / (f4t7 * rhydromod));

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			part = - hb_energy * f4t8Dsin / f4t8;
			Ftmp -= (a3 - rhydrodir * cosf(t8)) * (-part / rhydromod);
		  	Ttmp += _cross<number, number4>(rhydrodir, a3) * part;

			Ttmp += _cross<number, number4>(ppos_base, Ftmp);

			Ftmp.w = hb_energy;
			F += Ftmp;
		}
	}
	// END HYDROGEN BONDING

	// CROSS STACKING
	number4 rcstack = rhydro;
	number rcstackmodsqr = rhydromodsqr;
	if(SQR(CRST_RCLOW) < rcstackmodsqr && rcstackmodsqr < SQR(CRST_RCHIGH)) {
	  	number rcstackmod = sqrtf(rcstackmodsqr);
	  	number4 rcstackdir = rcstack / rcstackmod;

		// angles involved in the CSTCK interaction
		number t1 = CUDA_LRACOS (-CUDA_DOT(a1, b1));
		number t2 = CUDA_LRACOS (-CUDA_DOT(b1, rcstackdir));
		number t3 = CUDA_LRACOS ( CUDA_DOT(a1, rcstackdir));
		number t4 = CUDA_LRACOS ( CUDA_DOT(a3, b3));
		number t7 = CUDA_LRACOS (-CUDA_DOT(rcstackdir, b3));
		number t8 = CUDA_LRACOS ( CUDA_DOT(rcstackdir, a3));

	 	 // functions called at their relevant arguments
		number f2 = _f2(rcstackmod, CRST_F2);
		number f4t1 = _f4(t1, CRST_THETA1_T0, CRST_THETA1_TS, CRST_THETA1_TC, CRST_THETA1_A, CRST_THETA1_B);
		number f4t2 = _f4(t2, CRST_THETA2_T0, CRST_THETA2_TS, CRST_THETA2_TC, CRST_THETA2_A, CRST_THETA2_B);
		number f4t3 = _f4(t3, CRST_THETA3_T0, CRST_THETA3_TS, CRST_THETA3_TC, CRST_THETA3_A, CRST_THETA3_B);
		number f4t4 = _f4(t4, CRST_THETA4_T0, CRST_THETA4_TS, CRST_THETA4_TC, CRST_THETA4_A, CRST_THETA4_B) +
				_f4(PI - t4, CRST_THETA4_T0, CRST_THETA4_TS, CRST_THETA4_TC, CRST_THETA4_A, CRST_THETA4_B);
		number f4t7 = _f4(t7, CRST_THETA7_T0, CRST_THETA7_TS, CRST_THETA7_TC, CRST_THETA7_A, CRST_THETA7_B) +
				_f4(PI - t7, CRST_THETA7_T0, CRST_THETA7_TS, CRST_THETA7_TC, CRST_THETA7_A, CRST_THETA7_B);
		number f4t8 = _f4(t8, CRST_THETA8_T0, CRST_THETA8_TS, CRST_THETA8_TC, CRST_THETA8_A, CRST_THETA8_B) +
				_f4(PI - t8, CRST_THETA8_T0, CRST_THETA8_TS, CRST_THETA8_TC, CRST_THETA8_A, CRST_THETA8_B);

		number cstk_energy = f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

		if(cstk_energy < (number) 0) {
			// derivatives called at the relevant arguments
			number f2D = _f2D(rcstackmod, CRST_F2);
			number f4t1Dsin = -_f4Dsin(t1, CRST_THETA1_T0, CRST_THETA1_TS, CRST_THETA1_TC, CRST_THETA1_A, CRST_THETA1_B);
			number f4t2Dsin = -_f4Dsin(t2, CRST_THETA2_T0, CRST_THETA2_TS, CRST_THETA2_TC, CRST_THETA2_A, CRST_THETA2_B);
			number f4t3Dsin = _f4Dsin(t3, CRST_THETA3_T0, CRST_THETA3_TS, CRST_THETA3_TC, CRST_THETA3_A, CRST_THETA3_B);
			number f4t4Dsin = _f4Dsin(t4, CRST_THETA4_T0, CRST_THETA4_TS, CRST_THETA4_TC, CRST_THETA4_A, CRST_THETA4_B) -
					_f4Dsin(PI - t4, CRST_THETA4_T0, CRST_THETA4_TS, CRST_THETA4_TC, CRST_THETA4_A, CRST_THETA4_B);
			number f4t7Dsin = -_f4Dsin(t7, CRST_THETA7_T0, CRST_THETA7_TS, CRST_THETA7_TC, CRST_THETA7_A, CRST_THETA7_B) +
					_f4Dsin(PI - t7, CRST_THETA7_T0, CRST_THETA7_TS, CRST_THETA7_TC, CRST_THETA7_A, CRST_THETA7_B);
			number f4t8Dsin = _f4Dsin(t8, CRST_THETA8_T0, CRST_THETA8_TS, CRST_THETA8_TC, CRST_THETA8_A, CRST_THETA8_B) -
					_f4Dsin(PI - t8, CRST_THETA8_T0, CRST_THETA8_TS, CRST_THETA8_TC, CRST_THETA8_A, CRST_THETA8_B);

			// RADIAL PART
			Ftmp = rcstackdir * (cstk_energy * f2D / f2);

			// THETA1; t1 = LRACOS (-a1 * b1);
			Ttmp -= _cross<number, number4>(a1, b1) * (-cstk_energy * f4t1Dsin / f4t1);

			// TETA2; t2 = LRACOS (-b1 * rhydrodir);
			Ftmp -= (b1 + rcstackdir * cosf(t2)) * (cstk_energy * f4t2Dsin / (f4t2 * rcstackmod));

			// TETA3; t3 = LRACOS (a1 * rhydrodir);
			number part = -cstk_energy * f4t3Dsin / f4t3;
			Ftmp -= (a1 - rcstackdir * cosf(t3)) * (-part / rcstackmod);
			Ttmp += _cross<number, number4>(rcstackdir, a1) * part;

			// TETA4; t4 = LRACOS (a3 * b3);
			Ttmp -= _cross<number, number4>(a3, b3) * (-cstk_energy * f4t4Dsin / f4t4);

			// THETA7; t7 = LRACOS (-rcsrackir * b3);
			Ftmp -= (b3 + rcstackdir * cosf(t7)) * (cstk_energy * f4t7Dsin / (f4t7 * rcstackmod));

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			part = -cstk_energy * f4t8Dsin / f4t8;
			Ftmp -= (a3 - rcstackdir * cosf(t8)) * (-part / rcstackmod);
			Ttmp += _cross<number, number4>(rcstackdir, a3) * part;

			Ttmp += _cross<number, number4>(ppos_base, Ftmp);

			Ftmp.w = cstk_energy;
			F += Ftmp;
		}
	}

	// COAXIAL STACKING
	if (use_oxDNA2_coaxial_stacking){
		number4 rstack = r + qpos_stack - ppos_stack;
		number rstackmodsqr = CUDA_DOT(rstack, rstack);
		if(SQR(CXST_RCLOW) < rstackmodsqr && rstackmodsqr < SQR(CXST_RCHIGH)) {
			number rstackmod = sqrtf(rstackmodsqr);
			number4 rstackdir = rstack / rstackmod;

			// angles involved in the CXST interaction
			number t1 = CUDA_LRACOS (-CUDA_DOT(a1, b1));
			number t4 = CUDA_LRACOS ( CUDA_DOT(a3, b3));
			number t5 = CUDA_LRACOS ( CUDA_DOT(a3, rstackdir));
			number t6 = CUDA_LRACOS (-CUDA_DOT(b3, rstackdir));

			// functions called at their relevant arguments
			number f2 = _f2(rstackmod, CXST_F2);
			number f4t1 = _f4(t1, CXST_THETA1_T0_OXDNA2, CXST_THETA1_TS, CXST_THETA1_TC, CXST_THETA1_A, CXST_THETA1_B) +
				_f4_pure_harmonic(t1, CXST_THETA1_SA, CXST_THETA1_SB);
			number f4t4 = _f4(t4, CXST_THETA4_T0, CXST_THETA4_TS, CXST_THETA4_TC, CXST_THETA4_A, CXST_THETA4_B);
			number f4t5 = _f4(t5, CXST_THETA5_T0, CXST_THETA5_TS, CXST_THETA5_TC, CXST_THETA5_A, CXST_THETA5_B) +
					_f4(PI - t5, CXST_THETA5_T0, CXST_THETA5_TS, CXST_THETA5_TC, CXST_THETA5_A, CXST_THETA5_B);
			number f4t6 = _f4(t6, CXST_THETA6_T0, CXST_THETA6_TS, CXST_THETA6_TC, CXST_THETA6_A, CXST_THETA6_B) +
					_f4(PI - t6, CXST_THETA6_T0, CXST_THETA6_TS, CXST_THETA6_TC, CXST_THETA6_A, CXST_THETA6_B);

			number cxst_energy = f2 * f4t1 * f4t4 * f4t5 * f4t6;

			if(cxst_energy < (number) 0) {
				// derivatives called at the relevant arguments
				number f2D = _f2D(rstackmod, CXST_F2);
				number f4t1Dsin = -_f4Dsin(t1, CXST_THETA1_T0_OXDNA2, CXST_THETA1_TS, CXST_THETA1_TC, CXST_THETA1_A, CXST_THETA1_B) -
					_f4Dsin_pure_harmonic(t1, CXST_THETA1_SA, CXST_THETA1_SB);
				number f4t4Dsin =  _f4Dsin(t4, CXST_THETA4_T0, CXST_THETA4_TS, CXST_THETA4_TC, CXST_THETA4_A, CXST_THETA4_B);
				number f4t5Dsin =  _f4Dsin(t5, CXST_THETA5_T0, CXST_THETA5_TS, CXST_THETA5_TC, CXST_THETA5_A, CXST_THETA5_B) -
						_f4Dsin(PI - t5, CXST_THETA5_T0, CXST_THETA5_TS, CXST_THETA5_TC, CXST_THETA5_A, CXST_THETA5_B);
				number f4t6Dsin = -_f4Dsin(t6, CXST_THETA6_T0, CXST_THETA6_TS, CXST_THETA6_TC, CXST_THETA6_A, CXST_THETA6_B) +
						_f4Dsin(PI - t6, CXST_THETA6_T0, CXST_THETA6_TS, CXST_THETA6_TC, CXST_THETA6_A, CXST_THETA6_B);

				// RADIAL PART
				Ftmp = rstackdir * (cxst_energy * f2D / f2);

				// THETA1; t1 = LRACOS (-a1 * b1);
				Ttmp -= _cross<number, number4>(a1, b1) * (-cxst_energy * f4t1Dsin / f4t1);

				// TETA4; t4 = LRACOS (a3 * b3);
				Ttmp -= _cross<number, number4>(a3, b3) * (-cxst_energy * f4t4Dsin / f4t4);

				// THETA5; t5 = LRACOS ( a3 * rstackdir);
				number part = cxst_energy * f4t5Dsin / f4t5;
				Ftmp -= (a3 - rstackdir * cosf(t5)) / rstackmod * part;
				Ttmp -= _cross<number, number4>(rstackdir, a3) * part;

				// THETA6; t6 = LRACOS (-b3 * rstackdir);
				Ftmp -= (b3 + rstackdir * cosf(t6)) * (cxst_energy * f4t6Dsin / (f4t6 * rstackmod));

				Ttmp += _cross<number, number4>(ppos_stack, Ftmp);

				Ftmp.w = cxst_energy;
				F += Ftmp;
			}
		}
	}
	else {
		number4 rstack = r + qpos_stack - ppos_stack;
		number rstackmodsqr = CUDA_DOT(rstack, rstack);
		if(SQR(CXST_RCLOW) < rstackmodsqr && rstackmodsqr < SQR(CXST_RCHIGH)) {
			number rstackmod = sqrtf(rstackmodsqr);
			number4 rstackdir = rstack / rstackmod;

			// angles involved in the CXST interaction
			number t1 = CUDA_LRACOS (-CUDA_DOT(a1, b1));
			number t4 = CUDA_LRACOS ( CUDA_DOT(a3, b3));
			number t5 = CUDA_LRACOS ( CUDA_DOT(a3, rstackdir));
			number t6 = CUDA_LRACOS (-CUDA_DOT(b3, rstackdir));

			// This is the position the backbone would have with major-minor grooves the same width.
			// We need to do this to implement different major-minor groove widths because rback is
			// used as a reference point for things that have nothing to do with the actual backbone
			// position (in this case, the coaxial stacking interaction).
			number4 rbackboneref = r + POS_BACK * b1 - POS_BACK * a1;
			number rbackrefmod = _module<number, number4>(rbackboneref);
			number4 rbackbonerefdir = rbackboneref / rbackrefmod;
			number cosphi3 = CUDA_DOT(rstackdir, (_cross<number, number4>(rbackbonerefdir, a1)));

			// functions called at their relevant arguments
			number f2 = _f2(rstackmod, CXST_F2);
			number f4t1 = _f4(t1, CXST_THETA1_T0_OXDNA, CXST_THETA1_TS, CXST_THETA1_TC, CXST_THETA1_A, CXST_THETA1_B) +
					_f4(2 * PI - t1, CXST_THETA1_T0_OXDNA, CXST_THETA1_TS, CXST_THETA1_TC, CXST_THETA1_A, CXST_THETA1_B);
			number f4t4 = _f4(t4, CXST_THETA4_T0, CXST_THETA4_TS, CXST_THETA4_TC, CXST_THETA4_A, CXST_THETA4_B);
			number f4t5 = _f4(t5, CXST_THETA5_T0, CXST_THETA5_TS, CXST_THETA5_TC, CXST_THETA5_A, CXST_THETA5_B) +
					_f4(PI - t5, CXST_THETA5_T0, CXST_THETA5_TS, CXST_THETA5_TC, CXST_THETA5_A, CXST_THETA5_B);
			number f4t6 = _f4(t6, CXST_THETA6_T0, CXST_THETA6_TS, CXST_THETA6_TC, CXST_THETA6_A, CXST_THETA6_B) +
					_f4(PI - t6, CXST_THETA6_T0, CXST_THETA6_TS, CXST_THETA6_TC, CXST_THETA6_A, CXST_THETA6_B);
			number f5cosphi3 = _f5(cosphi3, CXST_F5_PHI3);

			number cxst_energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);

			if(cxst_energy < (number) 0) {
				// derivatives called at the relevant arguments
				number f2D = _f2D(rstackmod, CXST_F2);
				number f4t1Dsin = -_f4Dsin(t1, CXST_THETA1_T0_OXDNA, CXST_THETA1_TS, CXST_THETA1_TC, CXST_THETA1_A, CXST_THETA1_B) +
						_f4Dsin(2 * PI - t1, CXST_THETA1_T0_OXDNA, CXST_THETA1_TS, CXST_THETA1_TC, CXST_THETA1_A, CXST_THETA1_B);
				number f4t4Dsin =  _f4Dsin(t4, CXST_THETA4_T0, CXST_THETA4_TS, CXST_THETA4_TC, CXST_THETA4_A, CXST_THETA4_B);
				number f4t5Dsin =  _f4Dsin(t5, CXST_THETA5_T0, CXST_THETA5_TS, CXST_THETA5_TC, CXST_THETA5_A, CXST_THETA5_B) -
						_f4Dsin(PI - t5, CXST_THETA5_T0, CXST_THETA5_TS, CXST_THETA5_TC, CXST_THETA5_A, CXST_THETA5_B);
				number f4t6Dsin = -_f4Dsin(t6, CXST_THETA6_T0, CXST_THETA6_TS, CXST_THETA6_TC, CXST_THETA6_A, CXST_THETA6_B) +
						_f4Dsin(PI - t6, CXST_THETA6_T0, CXST_THETA6_TS, CXST_THETA6_TC, CXST_THETA6_A, CXST_THETA6_B);
				number f5cosphi3D = _f5D(cosphi3, CXST_F5_PHI3);

				// RADIAL PART
				Ftmp = rstackdir * (cxst_energy * f2D / f2);

				// THETA1; t1 = LRACOS (-a1 * b1);
				Ttmp -= _cross<number, number4>(a1, b1) * (-cxst_energy * f4t1Dsin / f4t1);

				// TETA4; t4 = LRACOS (a3 * b3);
				Ttmp -= _cross<number, number4>(a3, b3) * (-cxst_energy * f4t4Dsin / f4t4);

				// THETA5; t5 = LRACOS ( a3 * rstackdir);
				number part = cxst_energy * f4t5Dsin / f4t5;
				Ftmp -= (a3 - rstackdir * cosf(t5)) / rstackmod * part;
				Ttmp -= _cross<number, number4>(rstackdir, a3) * part;

				// THETA6; t6 = LRACOS (-b3 * rstackdir);
				Ftmp -= (b3 + rstackdir * cosf(t6)) * (cxst_energy * f4t6Dsin / (f4t6 * rstackmod));

				// COSPHI3
				number rbackrefmodcub = rbackrefmod * rbackrefmod * rbackrefmod;

				//number a1b1 = a1 * b1;
				number a2b1 = CUDA_DOT(a2, b1);
				number a3b1 = CUDA_DOT(a3, b1);
				number ra1 = CUDA_DOT(rstackdir, a1);
				number ra2 = CUDA_DOT(rstackdir, a2);
				number ra3 = CUDA_DOT(rstackdir, a3);
				number rb1 = CUDA_DOT(rstackdir, b1);

				number parentesi = (ra3 * a2b1 - ra2 * a3b1);
				number dcdr    = -GAMMA * parentesi * (GAMMA * (ra1 - rb1) + rstackmod) / rbackrefmodcub;
				number dcda1b1 =  GAMMA * SQR(GAMMA) * parentesi / rbackrefmodcub;
				number dcda2b1 =  GAMMA * ra3 / rbackrefmod;
				number dcda3b1 = -GAMMA * ra2 / rbackrefmod;
				number dcdra1  = -SQR(GAMMA) * parentesi * rstackmod / rbackrefmodcub;
				number dcdra2  = -GAMMA * a3b1 / rbackrefmod;
				number dcdra3  =  GAMMA * a2b1 / rbackrefmod;
				number dcdrb1  = -dcdra1;

				part = cxst_energy * 2 * f5cosphi3D / f5cosphi3;

				Ftmp -= part * (rstackdir * dcdr +
							    ((a1 - rstackdir * ra1) * dcdra1 +
								(a2 - rstackdir * ra2) * dcdra2 +
								(a3 - rstackdir * ra3) * dcdra3 +
								(b1 - rstackdir * rb1) * dcdrb1) / rstackmod);

				Ttmp += part * (_cross<number, number4>(rstackdir, a1) * dcdra1 +
							    _cross<number, number4>(rstackdir, a2) * dcdra2 +
							    _cross<number, number4>(rstackdir ,a3) * dcdra3);

				Ttmp -= part * (_cross<number, number4>(a1, b1) * dcda1b1 +
							    _cross<number, number4>(a2, b1) * dcda2b1 +
							    _cross<number, number4>(a3, b1) * dcda3b1);

				Ttmp += _cross<number, number4>(ppos_stack, Ftmp);

				Ftmp.w = cxst_energy;
				F += Ftmp;
			}
		}
	}
	
	// DEBYE HUCKEL
	if (use_debye_huckel){
		number rbackmod = _module<number, number4>(rbackbone);
		if (rbackmod < MD_dh_RC[0]){
			number4 rbackdir = rbackbone / rbackmod;
			if(rbackmod < MD_dh_RHIGH[0]){
				Ftmp = rbackdir * (-MD_dh_prefactor[0] * expf(MD_dh_minus_kappa[0] * rbackmod) * (MD_dh_minus_kappa[0] / rbackmod - 1.0f / SQR(rbackmod)));
			}
			else {
				Ftmp = rbackdir * (-2.0f * MD_dh_B[0] * (rbackmod - MD_dh_RC[0]));
			}

			// check for half-charge strand ends
			if (MD_dh_half_charged_ends[0] && (pbonds.n3 == P_INVALID || pbonds.n5 == P_INVALID)) {
				Ftmp *= 0.5f;
			}
			if (MD_dh_half_charged_ends[0] && (qbonds.n3 == P_INVALID || qbonds.n5 == P_INVALID)) {
				Ftmp *= 0.5f;
			}

			Ttmp -= _cross<number, number4>(ppos_back, Ftmp);
			F -= Ftmp;
		}
	}
	
	T += Ttmp;
	
	// this component stores the energy due to hydrogen bonding
	T.w = old_Tw + hb_energy;
}

// forces + second step without lists
template <typename number, typename number4>
__global__ void dna_forces(number4 *poss, GPU_quat<number> *orientations, number4 *forces, number4 *torques, LR_bonds *bonds, bool grooving, bool use_debye_huckel, bool use_oxDNA2_coaxial_stacking, bool use_oxDNA2_FENE) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	number4 T = make_number4<number, number4>(0, 0, 0, 0);
	LR_bonds bs = bonds[IND];
	number4 ppos = poss[IND];

	number4 a1, a2, a3;
	get_vectors_from_quat<number,number4>(orientations[IND], a1, a2, a3); //Returns vectors a1,a2 and a3 as they would be in the GPU matrix. These are necessary even in pure quaternion dynamics

	if(bs.n3 != P_INVALID) {
		number4 qpos = poss[bs.n3];
		number4 b1, b2, b3;
		get_vectors_from_quat<number,number4>(orientations[bs.n3], b1, b2, b3);

		_bonded_part<number, number4, true>(ppos, a1, a2, a3,
						    qpos, b1, b2, b3, F, T, grooving, use_oxDNA2_FENE);
	}

	if(bs.n5 != P_INVALID) {
		number4 qpos = poss[bs.n5];

		number4 b1, b2, b3;
		get_vectors_from_quat<number,number4>(orientations[bs.n5], b1, b2, b3);

		_bonded_part<number, number4, false>(qpos, b1, b2, b3,
						     ppos, a1, a2, a3, F, T, grooving, use_oxDNA2_FENE);
	}

	const int type = get_particle_type<number, number4>(ppos);
	T.w = (number) 0;
	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND && bs.n3 != j && bs.n5 != j) {
			const number4 qpos = poss[j];
			number4 b1, b2, b3;
			get_vectors_from_quat<number,number4>(orientations[j], b1, b2, b3);
			LR_bonds qbonds = bonds[j];

			_particle_particle_interaction<number, number4>(ppos, a1, a2, a3, qpos, b1, b2, b3, F, T, grooving, use_debye_huckel, use_oxDNA2_coaxial_stacking, bs, qbonds, IND, j);
		}
	}

	T = _vectors_transpose_number4_product(a1, a2, a3, T);

	// the real energy per particle is half of the one computed (because we count each interaction twice)
	F.w *= (number) 0.5f;
	T.w *= (number) 0.5f;
	forces[IND] = F;
	torques[IND] = T;
}

template <typename number, typename number4>
__global__ void dna_forces_edge_nonbonded(number4 *poss, GPU_quat<number> *orientations, number4 *forces, number4 *torques, edge_bond *edge_list, int n_edges, LR_bonds *bonds, bool grooving, bool use_debye_huckel, bool use_oxDNA2_coaxial_stacking) {
	if(IND >= n_edges) return;

	number4 dF = make_number4<number, number4>(0, 0, 0, 0);
	number4 dT = make_number4<number, number4>(0, 0, 0, 0);

	edge_bond b = edge_list[IND];

	// get info for particle 1
	number4 ppos = poss[b.from];
	// particle axes according to Allen's paper
	number4 a1, a2, a3;
	get_vectors_from_quat<number,number4>(orientations[b.from], a1, a2, a3);
	// get info for particle 2
	number4 qpos = poss[b.to];

	number4 b1, b2, b3;
	get_vectors_from_quat<number,number4>(orientations[b.to], b1, b2, b3);
	LR_bonds pbonds = bonds[b.from];
	LR_bonds qbonds = bonds[b.to];
	_particle_particle_interaction<number, number4>(ppos, a1, a2, a3, qpos, b1, b2, b3, dF, dT, grooving, use_debye_huckel, use_oxDNA2_coaxial_stacking, pbonds, qbonds, b.from, b.to);

	dF.w *= (number) 0.5f;
	dT.w *= (number) 0.5f;

	int from_index = MD_N[0]*(IND % MD_n_forces[0]) + b.from;
	//int from_index = MD_N[0]*(b.n_from % MD_n_forces[0]) + b.from;
	if((dF.x*dF.x + dF.y*dF.y + dF.z*dF.z + dF.w*dF.w) > (number)0.f) LR_atomicAddXYZ(&(forces[from_index]), dF);
	if((dT.x*dT.x + dT.y*dT.y + dT.z*dT.z + dT.w*dT.w) > (number)0.f) LR_atomicAddXYZ(&(torques[from_index]), dT);

	// Allen Eq. 6 pag 3:
	number4 dr = minimum_image<number, number4>(ppos, qpos); // returns qpos-ppos
	number4 crx = _cross<number, number4> (dr, dF);
	dT.x = -dT.x + crx.x;
	dT.y = -dT.y + crx.y;
	dT.z = -dT.z + crx.z;

	dF.x = -dF.x;
	dF.y = -dF.y;
	dF.z = -dF.z;

	int to_index = MD_N[0]*(IND % MD_n_forces[0]) + b.to;
	//int to_index = MD_N[0]*(b.n_to % MD_n_forces[0]) + b.to;
	if((dF.x*dF.x + dF.y*dF.y + dF.z*dF.z + dF.w*dF.w) > (number)0.f) LR_atomicAddXYZ(&(forces[to_index]), dF);
	if((dT.x*dT.x + dT.y*dT.y + dT.z*dT.z + dT.w*dT.w) > (number)0.f) LR_atomicAddXYZ(&(torques[to_index]), dT);
}

// bonded interactions for edge-based approach
template <typename number, typename number4>
__global__ void dna_forces_edge_bonded(number4 *poss, GPU_quat<number> *orientations,  number4 *forces, number4 *torques, LR_bonds *bonds, bool grooving, bool use_oxDNA2_FENE) {
	if(IND >= MD_N[0]) return;

	number4 F0, T0;

	F0.x = forces[IND].x;
	F0.y = forces[IND].y;
	F0.z = forces[IND].z;
	F0.w = forces[IND].w;
	T0.x = torques[IND].x;
	T0.y = torques[IND].y;
	T0.z = torques[IND].z;
	T0.w = torques[IND].w;

	number4 dF = make_number4<number, number4>(0, 0, 0, 0);
	number4 dT = make_number4<number, number4>(0, 0, 0, 0);
	number4 ppos = poss[IND];
	LR_bonds bs = bonds[IND];
	// particle axes according to Allen's paper

	number4 a1, a2, a3;
	get_vectors_from_quat<number,number4>(orientations[IND], a1, a2, a3);

	if(bs.n3 != P_INVALID) {
		number4 qpos = poss[bs.n3];

		number4 b1,b2,b3;
		get_vectors_from_quat<number,number4>(orientations[bs.n3], b1, b2, b3);

		_bonded_part<number, number4, true>(ppos, a1, a2, a3, qpos, b1, b2, b3, dF, dT, grooving, use_oxDNA2_FENE);
	}
	if(bs.n5 != P_INVALID) {
		number4 qpos = poss[bs.n5];

		number4 b1,b2,b3;
		get_vectors_from_quat<number,number4>(orientations[bs.n5], b1, b2, b3);
		_bonded_part<number, number4, false>(qpos, b1, b2, b3, ppos, a1, a2, a3, dF, dT, grooving, use_oxDNA2_FENE);
	}

	// the real energy per particle is half of the one computed (because we count each interaction twice)
	dF.w *= (number) 0.5f;
	dT.w *= (number) 0.5f;

	forces[IND] = (dF + F0);
	torques[IND] = (dT + T0);

	torques[IND] = _vectors_transpose_number4_product(a1, a2, a3, torques[IND]);
}

// forces + second step with verlet lists
template <typename number, typename number4>
__global__ void dna_forces(number4 *poss, GPU_quat<number> *orientations,  number4 *forces, number4 *torques, int *matrix_neighs, int *number_neighs, LR_bonds *bonds, bool grooving, bool use_debye_huckel, bool use_oxDNA2_coaxial_stacking, bool use_oxDNA2_FENE) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	number4 T = make_number4<number, number4>(0, 0, 0, 0);
	number4 ppos = poss[IND];
	LR_bonds bs = bonds[IND];

	// particle axes according to Allen's paper
	number4 a1, a2, a3;
	get_vectors_from_quat<number,number4>(orientations[IND], a1, a2, a3);
	
	if(bs.n3 != P_INVALID) {
		number4 qpos = poss[bs.n3];
		number4 b1, b2, b3;
		get_vectors_from_quat<number,number4>(orientations[bs.n3], b1, b2, b3);
		_bonded_part<number, number4, true>(ppos, a1, a2, a3,
						    qpos, b1, b2, b3, F, T, grooving, use_oxDNA2_FENE);
	}
	if(bs.n5 != P_INVALID) {
		number4 qpos = poss[bs.n5];
		number4 b1,b2,b3;
		get_vectors_from_quat<number,number4>(orientations[bs.n5], b1, b2, b3);
		_bonded_part<number, number4, false>(qpos, b1, b2, b3,
						     ppos, a1, a2, a3, F, T, grooving, use_oxDNA2_FENE);
	}

	const int type = get_particle_type<number, number4>(ppos);
	const int num_neighs = number_neighs[IND];

	T.w = (number) 0;
	for(int j = 0; j < num_neighs; j++) {
		const int k_index = matrix_neighs[j*MD_N[0] + IND];

		const number4 qpos = poss[k_index];
		number4 b1, b2, b3;
		get_vectors_from_quat<number,number4>(orientations[k_index], b1, b2, b3);
		LR_bonds pbonds = bonds[IND];
		LR_bonds qbonds = bonds[k_index];
		_particle_particle_interaction<number, number4>(ppos, a1, a2, a3, qpos, b1, b2, b3, F, T, grooving, use_debye_huckel, use_oxDNA2_coaxial_stacking, pbonds, qbonds, IND, k_index);
	}
	
	T = _vectors_transpose_number4_product(a1, a2, a3, T);

	// the real energy per particle is half of the one computed (because we count each interaction twice)
	F.w *= (number) 0.5f;
	T.w *= (number) 0.5f;
	forces[IND] = F;
	torques[IND] = T;
}


//FFS order parameter pre-calculations

// check whether a particular pair of particles have hydrogen bonding energy lower than a given threshold hb_threshold (which may vary)
template <typename number, typename number4>
__global__ void hb_op_precalc(number4 *poss, GPU_quat<number> *orientations, int *op_pairs1, int *op_pairs2, float *hb_energies, int n_threads, bool *region_is_nearhb)
{
	if(IND >= n_threads) return;
	
	int pind = op_pairs1[IND];
	int qind = op_pairs2[IND];
	// get distance between this nucleotide pair's "com"s
	number4 ppos = poss[pind];
	number4 qpos = poss[qind];
	number4 r = minimum_image<number, number4>(ppos, qpos);

	// check whether hb energy is below a certain threshold for this nucleotide pair
	int ptype = get_particle_type<number, number4>(ppos);
	int qtype = get_particle_type<number, number4>(qpos);
	int pbtype = get_particle_btype<number, number4>(ppos);
	int qbtype = get_particle_btype<number, number4>(qpos);
	int int_type = pbtype + qbtype;

	GPU_quat<number> po = orientations[pind];
	GPU_quat<number> qo = orientations[qind];

	//This gets an extra two vectors that are not needed, but the function doesn't seem to be called at all, so should make no difference. 
	number4 a1, a2, a3, b1, b2, b3; 
	get_vectors_from_quat<number,number4>(po, a1, a2, a3);

	get_vectors_from_quat<number,number4>(qo, b1, b2, b3);

	number4 ppos_base = POS_BASE * a1;
	number4 qpos_base = POS_BASE * b1;

	// HYDROGEN BONDING
	number hb_energy = (number) 0;
	number4 rhydro = r + qpos_base - ppos_base;
	number rhydromodsqr = CUDA_DOT(rhydro, rhydro);
	if(int_type == 3 && SQR(HYDR_RCLOW) < rhydromodsqr && rhydromodsqr < SQR(HYDR_RCHIGH)) {
		// versor and magnitude of the base-base separation
	  	number rhydromod = sqrtf(rhydromodsqr);
	  	number4 rhydrodir = rhydro / rhydromod;

		// angles involved in the HB interaction
		number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));
		number t2 = CUDA_LRACOS(-CUDA_DOT(b1, rhydrodir));
		number t3 = CUDA_LRACOS(CUDA_DOT(a1, rhydrodir));
		number t4 = CUDA_LRACOS(CUDA_DOT(a3, b3));
		number t7 = CUDA_LRACOS(-CUDA_DOT(rhydrodir, b3));
		number t8 = CUDA_LRACOS(CUDA_DOT(rhydrodir, a3));

	 	 // functions called at their relevant arguments
		number f1 = _f1(rhydromod, HYDR_F1, ptype, qtype);
		number f4t1 = _f4(t1, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC, HYDR_THETA1_A, HYDR_THETA1_B);
		number f4t2 = _f4(t2, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC, HYDR_THETA2_A, HYDR_THETA2_B);
		number f4t3 = _f4(t3, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC, HYDR_THETA3_A, HYDR_THETA3_B);
		number f4t4 = _f4(t4, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC, HYDR_THETA4_A, HYDR_THETA4_B);
		number f4t7 = _f4(t7, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC, HYDR_THETA7_A, HYDR_THETA7_B);
		number f4t8 = _f4(t8, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC, HYDR_THETA8_A, HYDR_THETA8_B);

		hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
	}
	// END HYDROGEN BONDING

	hb_energies[IND] = hb_energy;
}

template <typename number, typename number4>
__global__ void near_hb_op_precalc(number4 *poss, GPU_quat<number> *orientations, int *op_pairs1, int *op_pairs2, bool *nearly_bonded_array, int n_threads, bool *region_is_nearhb)
{
	if(IND >= n_threads) return;
	
	int pind = op_pairs1[IND];
	int qind = op_pairs2[IND];
	// get distance between this nucleotide pair's "com"s
	number4 ppos = poss[pind];
	number4 qpos = poss[qind];
	number4 r = minimum_image<number, number4>(ppos, qpos);

	// check whether hb energy is below a certain threshold for this nucleotide pair
	int ptype = get_particle_type<number, number4>(ppos);
	int qtype = get_particle_type<number, number4>(qpos);
	int pbtype = get_particle_btype<number, number4>(ppos);
	int qbtype = get_particle_btype<number, number4>(qpos);
	int int_type = pbtype + qbtype;

	GPU_quat<number> po = orientations[pind];
	GPU_quat<number> qo = orientations[qind];

	//This gets extra a2 and b2 vectors that aren't needed. get_vectors_from_quat could easily be modified to only return the relevant vectors, but it will make computationally very little difference, since most of the same numbers need to be calculated anyway. Perhaps worth changing for memory considerations, however
	number4 a1, a2, a3, b1, b2, b3; 
	get_vectors_from_quat<number,number4>(po, a1, a2, a3);
	get_vectors_from_quat<number,number4>(qo, b1, b2, b3);

	number4 ppos_base = POS_BASE * a1;
	number4 qpos_base = POS_BASE * b1;

	// HYDROGEN BONDING
	number4 rhydro = r + qpos_base - ppos_base;
	number rhydromodsqr = CUDA_DOT(rhydro, rhydro);
	
	int total_nonzero;
	bool nearly_bonded = false;
	if(int_type == 3 && SQR(HYDR_RCLOW) < rhydromodsqr && rhydromodsqr < SQR(HYDR_RCHIGH)) {
		// versor and magnitude of the base-base separation
	  	number rhydromod = sqrtf(rhydromodsqr);
	  	number4 rhydrodir = rhydro / rhydromod;

		// angles involved in the HB interaction
		number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));
		number t2 = CUDA_LRACOS(-CUDA_DOT(b1, rhydrodir));
		number t3 = CUDA_LRACOS(CUDA_DOT(a1, rhydrodir));
		number t4 = CUDA_LRACOS(CUDA_DOT(a3, b3));
		number t7 = CUDA_LRACOS(-CUDA_DOT(rhydrodir, b3));
		number t8 = CUDA_LRACOS(CUDA_DOT(rhydrodir, a3));

		// functions called at their relevant arguments
		number f1 = _f1(rhydromod, HYDR_F1, ptype, qtype);
		number f4t1 = _f4(t1, HYDR_THETA1_T0, HYDR_THETA1_TS, HYDR_THETA1_TC, HYDR_THETA1_A, HYDR_THETA1_B);
		number f4t2 = _f4(t2, HYDR_THETA2_T0, HYDR_THETA2_TS, HYDR_THETA2_TC, HYDR_THETA2_A, HYDR_THETA2_B);
		number f4t3 = _f4(t3, HYDR_THETA3_T0, HYDR_THETA3_TS, HYDR_THETA3_TC, HYDR_THETA3_A, HYDR_THETA3_B);
		number f4t4 = _f4(t4, HYDR_THETA4_T0, HYDR_THETA4_TS, HYDR_THETA4_TC, HYDR_THETA4_A, HYDR_THETA4_B);
		number f4t7 = _f4(t7, HYDR_THETA7_T0, HYDR_THETA7_TS, HYDR_THETA7_TC, HYDR_THETA7_A, HYDR_THETA7_B);
		number f4t8 = _f4(t8, HYDR_THETA8_T0, HYDR_THETA8_TS, HYDR_THETA8_TC, HYDR_THETA8_A, HYDR_THETA8_B);

		// the nearly bonded order parameter requires either all or all but one of these factors to be non-zero
		bool f1not0 = f1 < 0 ? 1 : 0;
		bool f4t1not0 = f4t1 > 0 ? 1 : 0;
		bool f4t2not0 = f4t2 > 0 ? 1 : 0;
		bool f4t3not0 = f4t3 > 0 ? 1 : 0;
		bool f4t4not0 = f4t4 > 0 ? 1 : 0;
		bool f4t7not0 = f4t7 > 0 ? 1 : 0;
		bool f4t8not0 = f4t8 > 0 ? 1 : 0;
		total_nonzero = f1not0 + f4t1not0 + f4t2not0 + f4t3not0 + f4t4not0 + f4t7not0 + f4t8not0;
		if (total_nonzero >= 6) nearly_bonded = true;
	}
	// END HYDROGEN BONDING

	nearly_bonded_array[IND] = nearly_bonded;
}

// compute the distance between a pair of particles
template <typename number, typename number4>
__global__ void dist_op_precalc(number4 *poss, GPU_quat<number> *orientations, int *op_pairs1, int *op_pairs2, number *op_dists, int n_threads)
{
	if(IND >= n_threads) return;

	int pind = op_pairs1[IND];
	int qind = op_pairs2[IND];

	// get distance between this nucleotide pair's "com"s
	number4 ppos = poss[pind];
	number4 qpos = poss[qind];
	number4 r = minimum_image<number, number4>(ppos, qpos);

	GPU_quat<number> po = orientations[pind];
	GPU_quat<number> qo = orientations[qind];

	//This gets extra a2 and b2 vectors that aren't needed. get_vectors_from_quat could easily be modified to only return the relevant vectors, but it will make computationally very little difference, since most of the same numbers need to be calculated anyway. Perhaps worth changing for memory considerations, however 
	number4 a1, a2, a3, b1, b2, b3; 
	get_vectors_from_quat<number,number4>(po, a1, a2, a3);
	get_vectors_from_quat<number,number4>(qo, b1, b2, b3);

	number4 ppos_base = POS_BASE * a1;
	number4 qpos_base = POS_BASE * b1;
	
	number4 rbase = r + qpos_base - ppos_base;
	op_dists[IND] = _module<number, number4>(rbase);
}
