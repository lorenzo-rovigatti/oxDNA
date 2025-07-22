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

// oxDNA3
__constant__ OxDNA3Params MD_fene_r0_SD[1];
__constant__ OxDNA3Params MD_fene_delta2_SD[1];
__constant__ OxDNA3Params MD_mbf_xmax_SD[1];
__constant__ float MD_fene_eps[1];

// global memory arrays. Probably we will have to find a better way to handle all these arrays (textures?).
__device__ OxDNA3Params MD_excl_s[7];
__device__ OxDNA3Params MD_excl_r[7];
__device__ OxDNA3Params MD_excl_b[7];
__device__ OxDNA3Params MD_excl_rc[7];

__device__ OxDNA3Params MD_F1_SD_A[2];
__device__ OxDNA3Params MD_F1_SD_RC[2];
__device__ OxDNA3Params MD_F1_SD_R0[2];
__device__ OxDNA3Params MD_F1_SD_BLOW[2];
__device__ OxDNA3Params MD_F1_SD_BHIGH[2];
__device__ OxDNA3Params MD_F1_SD_RLOW[2];
__device__ OxDNA3Params MD_F1_SD_RHIGH[2];
__device__ OxDNA3Params MD_F1_SD_RCLOW[2];
__device__ OxDNA3Params MD_F1_SD_RCHIGH[2];
__device__ OxDNA3Params MD_F1_SD_SHIFT[2];

__device__ OxDNA3Params MD_F2_SD_K[4];
__device__ OxDNA3Params MD_F2_SD_RC[4];
__device__ OxDNA3Params MD_F2_SD_R0[4];
__device__ OxDNA3Params MD_F2_SD_BLOW[4];
__device__ OxDNA3Params MD_F2_SD_RLOW[4];
__device__ OxDNA3Params MD_F2_SD_RCLOW[4];
__device__ OxDNA3Params MD_F2_SD_BHIGH[4];
__device__ OxDNA3Params MD_F2_SD_RCHIGH[4];
__device__ OxDNA3Params MD_F2_SD_RHIGH[4];

__device__ OxDNA3Params MD_F4_SD_THETA_A[21];
__device__ OxDNA3Params MD_F4_SD_THETA_B[21];
__device__ OxDNA3Params MD_F4_SD_THETA_T0[21];
__device__ OxDNA3Params MD_F4_SD_THETA_TS[21];
__device__ OxDNA3Params MD_F4_SD_THETA_TC[21];

__device__ OxDNA3Params MD_F5_SD_PHI_A[4];
__device__ OxDNA3Params MD_F5_SD_PHI_B[4];
__device__ OxDNA3Params MD_F5_SD_PHI_XC[4];
__device__ OxDNA3Params MD_F5_SD_PHI_XS[4];

#include "../cuda_utils/CUDA_lr_common.cuh"

__device__ struct neigh_types {
	uint8_t n3, n5;
	__device__ neigh_types(uint8_t nn3, uint8_t nn5) : n3(nn3), n5(nn5) {}
	__device__ neigh_types(uint8_t *particle_types, LR_bonds &bonds) {
		n3 = (bonds.n3 == P_INVALID) ? NO_TYPE : particle_types[bonds.n3];
		n5 = (bonds.n5 == P_INVALID) ? NO_TYPE : particle_types[bonds.n5];
	}

	__device__ bool is_end() {
		return n3 == NO_TYPE || n5 == NO_TYPE;
	}
};

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

__forceinline__ __device__ c_number _f1_SD(c_number r, int type, int n3_2, int n3_1, int n5_1, int n5_2) {
	c_number val = (c_number) 0.f;
	if(r < MD_F1_SD_RCHIGH[type](n3_2, n3_1, n5_1, n5_2)) {
		int eps_index = 25 * type + n3_1 * 5 + n5_1;
		if(r > MD_F1_SD_RHIGH[type](n3_2, n3_1, n5_1, n5_2)) {
			val = MD_F1_EPS[eps_index] * MD_F1_SD_BHIGH[type](n3_2, n3_1, n5_1, n5_2) * SQR(r - MD_F1_SD_RCHIGH[type](n3_2, n3_1, n5_1, n5_2));
		}
		else if(r > MD_F1_SD_RLOW[type](n3_2, n3_1, n5_1, n5_2)) {
			c_number tmp = 1.f - expf(-(r - MD_F1_SD_R0[type](n3_2, n3_1, n5_1, n5_2)) * MD_F1_SD_A[type](n3_2, n3_1, n5_1, n5_2));
			val = MD_F1_EPS[eps_index] * SQR(tmp) - MD_F1_SD_SHIFT[type](n3_2, n3_1, n5_1, n5_2);
		}
		else if(r > MD_F1_SD_RCLOW[type](n3_2, n3_1, n5_1, n5_2)) {
			val = MD_F1_EPS[eps_index] * MD_F1_SD_BLOW[type](n3_2, n3_1, n5_1, n5_2) * SQR(r - MD_F1_SD_RCLOW[type](n3_2, n3_1, n5_1, n5_2));
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

__forceinline__ __device__ c_number _f1D_SD(c_number r, int type, int n3_2, int n3_1, int n5_1, int n5_2) {
	c_number val = 0.f;
	if(r < MD_F1_SD_RCHIGH[type](n3_2, n3_1, n5_1, n5_2)) {
		float eps = MD_F1_EPS[25 * type + n3_1 * 5 + n5_1];
		if(r > MD_F1_SD_RHIGH[type](n3_2, n3_1, n5_1, n5_2)) {
			val = 2.f * eps * MD_F1_SD_BHIGH[type](n3_2, n3_1, n5_1, n5_2) * (r - MD_F1_SD_RCHIGH[type](n3_2, n3_1, n5_1, n5_2));
		}
		else if(r > MD_F1_SD_RLOW[type](n3_2, n3_1, n5_1, n5_2)) {
			c_number tmp = expf(-(r - MD_F1_SD_R0[type](n3_2, n3_1, n5_1, n5_2)) * MD_F1_SD_A[type](n3_2, n3_1, n5_1, n5_2));
			val = 2.f * eps * (1.f - tmp) * tmp * MD_F1_SD_A[type](n3_2, n3_1, n5_1, n5_2);
		}
		else if(r > MD_F1_SD_RCLOW[type](n3_2, n3_1, n5_1, n5_2)) {
			val = 2.f * eps * MD_F1_SD_BLOW[type](n3_2, n3_1, n5_1, n5_2) * (r - MD_F1_SD_RCLOW[type](n3_2, n3_1, n5_1, n5_2));
		}
	}

	return val;
}

__forceinline__ __device__ c_number _f2(c_number r, int type) {
	c_number val = 0.f;
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

__forceinline__ __device__ c_number _f2_SD(c_number r, int type, int n3_2, int n3_1, int n5_1, int n5_2) {
	c_number val = 0.f;
    if(r < MD_F2_SD_RCHIGH[type](n3_2, n3_1, n5_1, n5_2)) {
        if(r > MD_F2_SD_RHIGH[type](n3_2, n3_1, n5_1, n5_2)) {
            val = MD_F2_SD_K[type](n3_2, n3_1, n5_1, n5_2) * MD_F2_SD_BHIGH[type](n3_2, n3_1, n5_1, n5_2) * SQR(r - MD_F2_SD_RCHIGH[type](n3_2, n3_1, n5_1, n5_2));
        }
		else if(r > MD_F2_SD_RLOW[type](n3_2, n3_1, n5_1, n5_2)) {
            val = (MD_F2_SD_K[type](n3_2, n3_1, n5_1, n5_2) / 2.) * (SQR(r - MD_F2_SD_R0[type](n3_2, n3_1, n5_1, n5_2)) - SQR(MD_F2_SD_RC[type](n3_2, n3_1, n5_1, n5_2) - MD_F2_SD_R0[type](n3_2, n3_1, n5_1, n5_2)));
        }
		else if(r > MD_F2_SD_RCLOW[type](n3_2, n3_1, n5_1, n5_2)) {
            val = MD_F2_SD_K[type](n3_2, n3_1, n5_1, n5_2) * MD_F2_SD_BLOW[type](n3_2, n3_1, n5_1, n5_2) * SQR(r - MD_F2_SD_RCLOW[type](n3_2, n3_1, n5_1, n5_2));
        }
    }
    return val;
}

__forceinline__ __device__ c_number _f2D_SD(number r, int type, int n3_2, int n3_1, int n5_1, int n5_2) {
    c_number val = 0.f;
    if(r < MD_F2_SD_RCHIGH[type](n3_2, n3_1, n5_1, n5_2)) {
        if(r > MD_F2_SD_RHIGH[type](n3_2, n3_1, n5_1, n5_2)) {
            val = 2.f * MD_F2_SD_K[type](n3_2, n3_1, n5_1, n5_2) * MD_F2_SD_BHIGH[type](n3_2, n3_1, n5_1, n5_2) * (r - MD_F2_SD_RCHIGH[type](n3_2, n3_1, n5_1, n5_2));
        } 
		else if(r > MD_F2_SD_RLOW[type](n3_2, n3_1, n5_1, n5_2)) {
            val = MD_F2_SD_K[type](n3_2, n3_1, n5_1, n5_2) * (r - MD_F2_SD_R0[type](n3_2, n3_1, n5_1, n5_2));
        } 
		else if(r > MD_F2_SD_RCLOW[type](n3_2, n3_1, n5_1, n5_2)) {
            val = 2.f * MD_F2_SD_K[type](n3_2, n3_1, n5_1, n5_2) * MD_F2_SD_BLOW[type](n3_2, n3_1, n5_1, n5_2) * (r - MD_F2_SD_RCLOW[type](n3_2, n3_1, n5_1, n5_2));
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

// TODO: this can be optimised as the one in oxDNA2
__forceinline__ __device__ c_number _f4_SD(c_number t, int type, int n3_2, int n3_1, int n5_1, int n5_2) {
	c_number val = 0.f;
    t -= MD_F4_SD_THETA_T0[type](n3_2, n3_1, n5_1, n5_2);
    if(t < 0)
        t *= -1;

    if(t < MD_F4_SD_THETA_TC[type](n3_2, n3_1, n5_1, n5_2)) {
        if(t > MD_F4_SD_THETA_TS[type](n3_2, n3_1, n5_1, n5_2)) {
            // smoothing
            val = MD_F4_SD_THETA_B[type](n3_2, n3_1, n5_1, n5_2) * SQR(MD_F4_SD_THETA_TC[type](n3_2, n3_1, n5_1, n5_2) - t);
        }
		else {
            val = 1.f - MD_F4_SD_THETA_A[type](n3_2, n3_1, n5_1, n5_2) * SQR(t);
		}
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

// TODO: this can be optimised as the one in oxDNA2
__forceinline__ __device__ c_number _f4D_SD(c_number t, int type, int n3_2, int n3_1, int n5_1, int n5_2) {
	c_number val = 0.f;
    c_number m = 1.f;
    t -= MD_F4_SD_THETA_T0[type](n3_2, n3_1, n5_1, n5_2);
    // this function is a parabola centered in t0. If t < 0 then the value of the function
    // is the same but the value of its derivative has the opposite sign, so m = -1
    if(t < 0) {
        t *= -1;
        m = -1.f;
    }

    if(t < MD_F4_SD_THETA_TC[type](n3_2, n3_1, n5_1, n5_2)) {
        if(t > MD_F4_SD_THETA_TS[type](n3_2, n3_1, n5_1, n5_2)) {
            // smoothing
            val = m * 2.f * MD_F4_SD_THETA_B[type](n3_2, n3_1, n5_1, n5_2) * (t - MD_F4_SD_THETA_TC[type](n3_2, n3_1, n5_1, n5_2));
        } else
            val = -m * 2.f * MD_F4_SD_THETA_A[type](n3_2, n3_1, n5_1, n5_2) * t;
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

__forceinline__ __device__ c_number _f5_SD(c_number f, int type, int n3_2, int n3_1, int n5_1, int n5_2) {
	c_number val = (c_number) 0.f;

	if(f > MD_F5_SD_PHI_XC[type](n3_2, n3_1, n5_1, n5_2)) {
		if(f < MD_F5_SD_PHI_XS[type](n3_2, n3_1, n5_1, n5_2)) {
			val = MD_F5_SD_PHI_B[type](n3_2, n3_1, n5_1, n5_2) * SQR(MD_F5_SD_PHI_XC[type](n3_2, n3_1, n5_1, n5_2) - f);
		}
		else if(f < 0.f) {
			val = (c_number) 1.f - MD_F5_SD_PHI_A[type](n3_2, n3_1, n5_1, n5_2) * SQR(f);
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

__forceinline__ __device__ c_number _f5D_SD(c_number f, int type, int n3_2, int n3_1, int n5_1, int n5_2) {
	c_number val = (c_number) 0.f;

	if(f > MD_F5_SD_PHI_XC[type](n3_2, n3_1, n5_1, n5_2)) {
		if(f < MD_F5_SD_PHI_XS[type](n3_2, n3_1, n5_1, n5_2)) {
			val = 2.f * MD_F5_SD_PHI_B[type](n3_2, n3_1, n5_1, n5_2) * (f - MD_F5_SD_PHI_XC[type](n3_2, n3_1, n5_1, n5_2));
		}
		else if(f < 0.f) {
			val = (c_number) -2.f * MD_F5_SD_PHI_A[type](n3_2, n3_1, n5_1, n5_2) * f;
		}
	}

	return val;
}

template<bool qIsN3>
__device__ void _DNA3_bonded_excluded_volume(const c_number4 &r, const c_number4 &n3pos_base, const c_number4 &n3pos_back, 
		int n3_type, int neigh_n3_type, const c_number4 &n5pos_base, const c_number4 &n5pos_back, int n5_type, int neigh_n5_type,
		c_number4 &F, c_number4 &T) {
	c_number4 Ftmp;
	// BASE-BASE
	c_number4 rcenter = r + n3pos_base - n5pos_base;
	c_number s = MD_excl_s[4](neigh_n3_type, n3_type, n5_type, neigh_n5_type);
	c_number rstar = MD_excl_r[4](neigh_n3_type, n3_type, n5_type, neigh_n5_type);
	c_number b = MD_excl_b[4](neigh_n3_type, n3_type, n5_type, neigh_n5_type);
	c_number rc = MD_excl_rc[4](neigh_n3_type, n3_type, n5_type, neigh_n5_type);
	_excluded_volume(rcenter, Ftmp, s, rstar, b, rc);
	c_number4 torquep1 = (qIsN3) ? _cross(n5pos_base, Ftmp) : _cross(n3pos_base, Ftmp);
	F += Ftmp;

	// n5-BASE vs. n3-BACK
	rcenter = r + n3pos_back - n5pos_base;
	s = MD_excl_s[5](neigh_n3_type, n3_type, n5_type, neigh_n5_type);
	rstar = MD_excl_r[5](neigh_n3_type, n3_type, n5_type, neigh_n5_type);
	b = MD_excl_b[5](neigh_n3_type, n3_type, n5_type, neigh_n5_type);
	rc = MD_excl_rc[5](neigh_n3_type, n3_type, n5_type, neigh_n5_type);
	_excluded_volume(rcenter, Ftmp, s, rstar, b, rc);
	c_number4 torquep2 = (qIsN3) ? _cross(n5pos_base, Ftmp) : _cross(n3pos_back, Ftmp);
	F += Ftmp;

	// n5-BACK vs. n3-BASE
	rcenter = r + n3pos_base - n5pos_back;
	s = MD_excl_s[6](neigh_n3_type, n3_type, n5_type, neigh_n5_type);
	rstar = MD_excl_r[6](neigh_n3_type, n3_type, n5_type, neigh_n5_type);
	b = MD_excl_b[6](neigh_n3_type, n3_type, n5_type, neigh_n5_type);
	rc = MD_excl_rc[6](neigh_n3_type, n3_type, n5_type, neigh_n5_type);
	_excluded_volume(rcenter, Ftmp, s, rstar, b, rc);
	c_number4 torquep3 = (qIsN3) ? _cross(n5pos_back, Ftmp) : _cross(n3pos_base, Ftmp);
	F += Ftmp;

	T += torquep1 + torquep2 + torquep3;
}

__device__ void _DNA3_nonbonded_excluded_volume(const c_number4 &r, const c_number4 &qpos_base, const c_number4 &qpos_back, 
		int qtype, int neigh_qtype, const c_number4 &ppos_base, const c_number4 &ppos_back, int ptype, int neigh_ptype,
		c_number4 &F, c_number4 &T) {
	c_number4 Ftmp;
	// BASE-BASE
	c_number4 rcenter = r + qpos_base - ppos_base;
	c_number s = MD_excl_s[1](neigh_qtype, qtype, ptype, neigh_ptype);
	c_number rstar = MD_excl_r[1](neigh_qtype, qtype, ptype, neigh_ptype);
	c_number b = MD_excl_b[1](neigh_qtype, qtype, ptype, neigh_ptype);
	c_number rc = MD_excl_rc[1](neigh_qtype, qtype, ptype, neigh_ptype);
	_excluded_volume(rcenter, Ftmp, s, rstar, b, rc);
	c_number4 torquep1 = _cross(ppos_base, Ftmp);
	F += Ftmp;

	// n5-BASE vs. n3-BACK
	rcenter = r + qpos_back - ppos_base;
	s = MD_excl_s[3](neigh_qtype, qtype, ptype, neigh_ptype);
	rstar = MD_excl_r[3](neigh_qtype, qtype, ptype, neigh_ptype);
	b = MD_excl_b[3](neigh_qtype, qtype, ptype, neigh_ptype);
	rc = MD_excl_rc[3](neigh_qtype, qtype, ptype, neigh_ptype);
	_excluded_volume(rcenter, Ftmp, s, rstar, b, rc);
	c_number4 torquep2 = _cross(ppos_base, Ftmp);
	F += Ftmp;

	// n5-BACK vs. n3-BASE
	rcenter = r + qpos_base - ppos_back;
	s = MD_excl_s[2](neigh_qtype, qtype, ptype, neigh_ptype);
	rstar = MD_excl_r[2](neigh_qtype, qtype, ptype, neigh_ptype);
	b = MD_excl_b[2](neigh_qtype, qtype, ptype, neigh_ptype);
	rc = MD_excl_rc[2](neigh_qtype, qtype, ptype, neigh_ptype);
	_excluded_volume(rcenter, Ftmp, s, rstar, b, rc);
	c_number4 torquep3 = _cross(ppos_back, Ftmp);
	F += Ftmp;

	// BACK-BACK
	rcenter = r + qpos_back - ppos_back;
	s = MD_excl_s[0](neigh_qtype, qtype, ptype, neigh_ptype);
	rstar = MD_excl_r[0](neigh_qtype, qtype, ptype, neigh_ptype);
	b = MD_excl_b[0](neigh_qtype, qtype, ptype, neigh_ptype);
	rc = MD_excl_rc[0](neigh_qtype, qtype, ptype, neigh_ptype);
	_excluded_volume(rcenter, Ftmp, s, rstar, b, rc);
	c_number4 torquep4 = _cross(ppos_back, Ftmp);
	F += Ftmp;

	T += torquep1 + torquep2 + torquep3 + torquep4;
}

__device__ void DNA3_set_interaction_sites(int type, int btype, const c_number4 &principal_axis, const c_number4 &third_axis, 
	c_number4 &r_back, c_number4 &r_stack, c_number4 &r_base) {
	c_number pos_mm_back1, pos_mm_back2, pos_stack, pos_base;

	if(btype == 4) {
		pos_mm_back1 = POS_MM_BACK1;
		pos_mm_back2 = POS_MM_BACK2;
		pos_stack = POS_STACK;
		pos_base =  POS_BASE;
	}
	else {
		if(type == N_A) {
			pos_mm_back1 = POS_MM_BACK1_A;
			pos_mm_back2 = POS_MM_BACK2_A;
			pos_stack = POS_STACK_A;
			pos_base =  POS_BASE_A;
		}
		else if(type == N_G) {
			pos_mm_back1 = POS_MM_BACK1_G;
			pos_mm_back2 = POS_MM_BACK2_G;
			pos_stack = POS_STACK_G;
			pos_base =  POS_BASE_G;
		}
		else if(type == N_C) {
			pos_mm_back1 = POS_MM_BACK1_C;
			pos_mm_back2 = POS_MM_BACK2_C;
			pos_stack = POS_STACK_C;
			pos_base =  POS_BASE_C;
		}
		else if(type == N_T) {
			pos_mm_back1 = POS_MM_BACK1_T;
			pos_mm_back2 = POS_MM_BACK2_T;
			pos_stack = POS_STACK_T;
			pos_base =  POS_BASE_T;
		}
	}
	
	r_back = principal_axis * pos_mm_back1 + third_axis * pos_mm_back2;
	r_stack = principal_axis * pos_stack;
	r_base = principal_axis * pos_base;
}

template<bool qIsN3>
__device__ void _DNA3_bonded_part(const c_number4 &r, const c_number4 &n5pos, const c_number4 &n5x, const c_number4 &n5y, const c_number4 &n5z, uint8_t neigh_n5_type, 
	const c_number4 &n3pos, const c_number4 &n3x, const c_number4 &n3y, const c_number4 &n3z, uint8_t neigh_n3_type, 
	c_number4 &F, c_number4 &T, bool use_mbf, c_number mbf_finf) {
	int n3type = get_particle_type(n3pos);
	int n5type = get_particle_type(n5pos);
	int n3btype = get_particle_btype(n3pos);
	int n5btype = get_particle_btype(n5pos);
	int n3_id = get_particle_index(n3pos);
	int n5_id = get_particle_index(n5pos);

	c_number4 n5pos_back, n5pos_base, n5pos_stack;
	DNA3_set_interaction_sites(n5type, n5btype, n5x, n5y, n5pos_back, n5pos_stack, n5pos_base);
	c_number4 n3pos_back, n3pos_base, n3pos_stack;
	DNA3_set_interaction_sites(n3type, n3btype, n3x, n3y, n3pos_back, n3pos_stack, n3pos_base);

	c_number4 rback = r + n3pos_back - n5pos_back;
	c_number rbackmod = _module(rback);
	c_number rbackr0 = rbackmod - MD_fene_r0_SD[0](neigh_n3_type, n3type, n5type, neigh_n5_type);

	c_number4 Ftmp;
	c_number mbf_xmax = MD_mbf_xmax_SD[0](neigh_n3_type, n3type, n5type, neigh_n5_type);
	c_number fene_delta2 = MD_fene_delta2_SD[0](neigh_n3_type, n3type, n5type, neigh_n5_type);
	if(use_mbf == true && fabsf(rbackr0) > mbf_xmax) {
		// this is the "relax" potential, i.e. the standard FENE up to xmax and then something like A + B log(r) for r>xmax
		c_number fene_xmax = -(MD_fene_eps[0] / 2.f) * logf(1.f - mbf_xmax * mbf_xmax / fene_delta2);
		c_number mbf_fmax = (MD_fene_eps[0] * mbf_xmax / (fene_delta2 - SQR(mbf_xmax)));
		c_number long_xmax = (mbf_fmax - mbf_finf) * mbf_xmax * logf(mbf_xmax) + mbf_finf * mbf_xmax;
		Ftmp = rback * (copysignf(1.f, rbackr0) * ((mbf_fmax - mbf_finf) * mbf_xmax / fabsf(rbackr0) + mbf_finf) / rbackmod);
		Ftmp.w = (mbf_fmax - mbf_finf) * mbf_xmax * logf(fabsf(rbackr0)) + mbf_finf * fabsf(rbackr0) - long_xmax + fene_xmax;
	}
	else {
		Ftmp = rback * ((MD_fene_eps[0] * rbackr0 / (fene_delta2 - SQR(rbackr0))) / rbackmod);
		Ftmp.w = -MD_fene_eps[0] * ((c_number) 0.5f) * logf(1 - SQR(rbackr0) / fene_delta2);
	}

	c_number4 Ttmp = (qIsN3) ? _cross(n5pos_back, Ftmp) : _cross(n3pos_back, Ftmp);
	// EXCLUDED VOLUME
	_DNA3_bonded_excluded_volume<qIsN3>(r, n3pos_base, n3pos_back, n3type, neigh_n3_type, n5pos_base, n5pos_back, n5type, neigh_n5_type, Ftmp, Ttmp);

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
	c_number f1 = _f1_SD(rstackmod, STCK_F1, neigh_n3_type, n3type, n5type, neigh_n5_type);
    c_number f4t4 = _f4_SD(t4, STCK_F4_THETA4, neigh_n3_type, n3type, n5type, neigh_n5_type);
    c_number f4t5 = _f4_SD(PI - t5, STCK_F4_THETA5, neigh_n3_type, n3type, n5type, neigh_n5_type);
    c_number f4t6 = _f4_SD(t6, STCK_F4_THETA6, neigh_n3_type, n3type, n5type, neigh_n5_type);
    c_number f5phi1 = _f5_SD(cosphi1, STCK_F5_PHI1, neigh_n3_type, n3type, n5type, neigh_n5_type);
    c_number f5phi2 = _f5_SD(cosphi2, STCK_F5_PHI2, neigh_n3_type, n3type, n5type, neigh_n5_type);

	c_number energy = f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;

	if(energy != (c_number) 0) {
		// and their derivatives
		c_number f1D = _f1D_SD(rstackmod, STCK_F1, neigh_n3_type, n3type, n5type, neigh_n5_type);
        c_number f4t4D = _f4D_SD(t4, STCK_F4_THETA4, neigh_n3_type, n3type, n5type, neigh_n5_type);
        c_number f4t5D = _f4D_SD(PI - t5, STCK_F4_THETA5, neigh_n3_type, n3type, n5type, neigh_n5_type);
        c_number f4t6D = _f4D_SD(t6, STCK_F4_THETA6, neigh_n3_type, n3type, n5type, neigh_n5_type);
        c_number f5phi1D = _f5D_SD(cosphi1, STCK_F5_PHI1, neigh_n3_type, n3type, n5type, neigh_n5_type);
        c_number f5phi2D = _f5D_SD(cosphi2, STCK_F5_PHI2, neigh_n3_type, n3type, n5type, neigh_n5_type);

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

__device__ void _DNA3_particle_particle_DNA_interaction(const c_number4 &r, const c_number4 &ppos, const c_number4 &a1, const c_number4 &a2, const c_number4 &a3,
		const c_number4 &qpos, const c_number4 &b1,	const c_number4 &b2, const c_number4 &b3, c_number4 &F, c_number4 &T, neigh_types &p_neighs, neigh_types &q_neighs){
	int ptype = get_particle_type(ppos);
	int qtype = get_particle_type(qpos);
	int pbtype = get_particle_btype(ppos);
	int qbtype = get_particle_btype(qpos);
	int int_type = pbtype + qbtype;

	c_number4 ppos_back, ppos_base, ppos_stack;
	DNA3_set_interaction_sites(ptype, pbtype, a1, a2, ppos_back, ppos_stack, ppos_base);
	c_number4 qpos_back, qpos_base, qpos_stack;
	DNA3_set_interaction_sites(qtype, qbtype, b1, b2, qpos_back, qpos_stack, qpos_base);
	c_number tot_Tw = T.w;

	// excluded volume
	c_number4 Ftmp = make_c_number4(0, 0, 0, 0);
	c_number4 Ttmp = make_c_number4(0, 0, 0, 0);
	_DNA3_nonbonded_excluded_volume(r, qpos_base, qpos_back, qtype, NO_TYPE, ppos_base, ppos_back, ptype, NO_TYPE, Ftmp, Ttmp);
	F += Ftmp;

	// DEBYE HUCKEL
	c_number4 rbackbone = r + qpos_back - ppos_back;
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
		if(MD_dh_half_charged_ends[0] && p_neighs.is_end()) {
			Ftmp *= 0.5f;
		}

		if(MD_dh_half_charged_ends[0] && q_neighs.is_end()) {
			Ftmp *= 0.5f;
		}

		Ttmp -= _cross(ppos_back, Ftmp);
		F -= Ftmp;
	}

	// HYDROGEN BONDING
	c_number4 rhydro = r + qpos_base - ppos_base;
	c_number rhydromodsqr = CUDA_DOT(rhydro, rhydro);
	if(int_type == 3 && SQR(MD_F1_SD_RCLOW[HYDR_F1](NO_TYPE, qtype, ptype, NO_TYPE)) < rhydromodsqr && rhydromodsqr < SQR(MD_F1_SD_RCHIGH[HYDR_F1](NO_TYPE, qtype, ptype, NO_TYPE))) {
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
		c_number f1 = hb_multi * _f1_SD(rhydromod, HYDR_F1, 0, qtype, ptype, 0);
		c_number f4t1 = _f4_SD(t1, HYDR_F4_THETA1, 0, qtype, ptype, 0);
		c_number f4t2 = _f4_SD(t2, HYDR_F4_THETA2, 0, qtype, ptype, 0);
		c_number f4t3 = _f4_SD(t3, HYDR_F4_THETA3, 0, qtype, ptype, 0);
		c_number f4t4 = _f4_SD(t4, HYDR_F4_THETA4, 0, qtype, ptype, 0);
		c_number f4t7 = _f4_SD(t7, HYDR_F4_THETA7, 0, qtype, ptype, 0);
		c_number f4t8 = _f4_SD(t8, HYDR_F4_THETA8, 0, qtype, ptype, 0);

		c_number hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
		tot_Tw += hb_energy;

		if(hb_energy < (c_number) 0) {
			// derivatives called at the relevant arguments
			c_number f1D = hb_multi * _f1D_SD(rhydromod, HYDR_F1, 0, qtype, ptype, 0);
			c_number f4t1D = -_f4D_SD(t1, HYDR_F4_THETA1, 0, qtype, ptype, 0);
			c_number f4t2D = -_f4D_SD(t2, HYDR_F4_THETA2, 0, qtype, ptype, 0);
			c_number f4t3D = _f4D_SD(t3, HYDR_F4_THETA3, 0, qtype, ptype, 0);
			c_number f4t4D = _f4D_SD(t4, HYDR_F4_THETA4, 0, qtype, ptype, 0);
			c_number f4t7D = -_f4D_SD(t7, HYDR_F4_THETA7, 0, qtype, ptype, 0);
			c_number f4t8D = _f4D_SD(t8, HYDR_F4_THETA8, 0, qtype, ptype, 0);

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
	int type_p_n3 = p_neighs.n3;
	int type_p_n5 = p_neighs.n5;
	int type_q_n3 = q_neighs.n3;
	int type_q_n5 = q_neighs.n5;
	c_number4 rcstack = rhydro;
	c_number rcstackmodsqr = rhydromodsqr;
	c_number rcstackmod = sqrtf(rcstackmodsqr);
	c_number4 rcstackdir = rcstack / rcstackmod;
	c_number cost7 = -CUDA_DOT(rcstackdir, b3);
	c_number cost8 = CUDA_DOT(rcstackdir, a3);
	if((cost7 > 0 && cost8 > 0 && MD_F2_SD_RCLOW[CRST_F2_33](type_q_n3, qtype, ptype, type_p_n3) < rcstackmod && rcstackmod < MD_F2_SD_RCHIGH[CRST_F2_33](type_q_n3, qtype, ptype, type_p_n3)) 
	|| (cost7 < 0 && cost8 < 0 && MD_F2_SD_RCLOW[CRST_F2_55](type_q_n5, qtype, ptype, type_p_n5) < rcstackmod && rcstackmod < MD_F2_SD_RCHIGH[CRST_F2_55](type_q_n5, qtype, ptype, type_p_n5))) {
		// angles involved in the CSTCK interaction
		c_number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));
		c_number cost2 = -CUDA_DOT(b1, rcstackdir);
		c_number t2 = CUDA_LRACOS(cost2);
		c_number cost3 = CUDA_DOT(a1, rcstackdir);
		c_number t3 = CUDA_LRACOS(cost3);
		c_number t4 = CUDA_LRACOS(CUDA_DOT(a3, b3));
		c_number t7 = CUDA_LRACOS(cost7);
		c_number t8 = CUDA_LRACOS(cost8);

		// functions called at their relevant arguments
		c_number f2_33 = _f2_SD(rcstackmod, CRST_F2_33, type_q_n3, qtype, ptype, type_p_n3);
        c_number f4t1_33 = _f4_SD(t1, CRST_F4_THETA1_33, type_q_n3, qtype, ptype, type_p_n3);
        c_number f4t2_33 = _f4_SD(t2, CRST_F4_THETA2_33, type_q_n3, qtype, ptype, type_p_n3);
        c_number f4t3_33 = _f4_SD(t3, CRST_F4_THETA3_33, type_q_n3, qtype, ptype, type_p_n3);
        c_number f4t4_33 = _f4_SD(t4, CRST_F4_THETA4_33, type_q_n3, qtype, ptype, type_p_n3);
        c_number f4t7_33 = _f4_SD(t7, CRST_F4_THETA7_33, type_q_n3, qtype, ptype, type_p_n3);
        c_number f4t8_33 = _f4_SD(t8, CRST_F4_THETA8_33, type_q_n3, qtype, ptype, type_p_n3);

        // 5'5' diagonal
        c_number f2_55 = _f2_SD(rcstackmod, CRST_F2_55, type_q_n5, qtype, ptype, type_p_n5);
        c_number f4t1_55 = _f4_SD(t1, CRST_F4_THETA1_55, type_q_n5, qtype, ptype, type_p_n5);
        c_number f4t2_55 = _f4_SD(t2, CRST_F4_THETA2_55, type_q_n5, qtype, ptype, type_p_n5);
        c_number f4t3_55 = _f4_SD(t3, CRST_F4_THETA3_55, type_q_n5, qtype, ptype, type_p_n5);
        c_number f4t4_55 = _f4_SD(t4, CRST_F4_THETA4_55, type_q_n5, qtype, ptype, type_p_n5);
        c_number f4t7_55 = _f4_SD(t7, CRST_F4_THETA7_55, type_q_n5, qtype, ptype, type_p_n5);
        c_number f4t8_55 = _f4_SD(t8, CRST_F4_THETA8_55, type_q_n5, qtype, ptype, type_p_n5);

		c_number cstk_energy = f2_33 * f4t1_33 * f4t2_33 * f4t3_33 * f4t4_33 * f4t7_33 * f4t8_33 + f2_55 * f4t1_55 * f4t2_55 * f4t3_55 * f4t4_55 * f4t7_55 * f4t8_55;

		if(cstk_energy < (c_number) 0) {
			// derivatives called at the relevant arguments
            number f2D_33 = _f2D_SD(rcstackmod, CRST_F2_33, type_q_n3, qtype, ptype, type_p_n3);
            number f4t1Dsin_33 = -_f4D_SD(t1, CRST_F4_THETA1_33, type_q_n3, qtype, ptype, type_p_n3);
            number f4t2Dsin_33 = -_f4D_SD(t2, CRST_F4_THETA2_33, type_q_n3, qtype, ptype, type_p_n3);
            number f4t3Dsin_33 =  _f4D_SD(t3, CRST_F4_THETA3_33, type_q_n3, qtype, ptype, type_p_n3);
            number f4t4Dsin_33 =  _f4D_SD(t4, CRST_F4_THETA4_33, type_q_n3, qtype, ptype, type_p_n3);
            number f4t7Dsin_33 = -_f4D_SD(t7, CRST_F4_THETA7_33, type_q_n3, qtype, ptype, type_p_n3);
            number f4t8Dsin_33 =  _f4D_SD(t8, CRST_F4_THETA8_33, type_q_n3, qtype, ptype, type_p_n3);

            number f2D_55 = _f2D_SD(rcstackmod, CRST_F2_55, type_q_n5, qtype, ptype, type_p_n5);
            number f4t1Dsin_55 = -_f4D_SD(t1, CRST_F4_THETA1_55, type_q_n5, qtype, ptype, type_p_n5);
            number f4t2Dsin_55 = -_f4D_SD(t2, CRST_F4_THETA2_55, type_q_n5, qtype, ptype, type_p_n5);
            number f4t3Dsin_55 =  _f4D_SD(t3, CRST_F4_THETA3_55, type_q_n5, qtype, ptype, type_p_n5);
            number f4t4Dsin_55 =  _f4D_SD(t4, CRST_F4_THETA4_55, type_q_n5, qtype, ptype, type_p_n5);
            number f4t7Dsin_55 = -_f4D_SD(t7, CRST_F4_THETA7_55, type_q_n5, qtype, ptype, type_p_n5);
            number f4t8Dsin_55 =  _f4D_SD(t8, CRST_F4_THETA8_55, type_q_n5, qtype, ptype, type_p_n5);

			// RADIAL PART
			Ftmp = rcstackdir * ((f2D_33 * f4t1_33 * f4t2_33 * f4t3_33 * f4t4_33 * f4t7_33 * f4t8_33) + (f2D_55 * f4t1_55 * f4t2_55 * f4t3_55 * f4t4_55 * f4t7_55 * f4t8_55));

			// THETA1; t1 = LRACOS (-a1 * b1);
			Ttmp -= stably_normalised(_cross(a1, b1)) * (-f2_33 * f4t1Dsin_33 * f4t2_33 * f4t3_33 * f4t4_33 * f4t7_33 * f4t8_33 - f2_55 * f4t1Dsin_55 * f4t2_55 * f4t3_55 * f4t4_55 * f4t7_55 * f4t8_55);

			// TETA2; t2 = LRACOS (-b1 * rhydrodir);
			Ftmp -= stably_normalised(b1 + rcstackdir * cost2) * ((f2_33 * f4t1_33 * f4t2Dsin_33 * f4t3_33 * f4t4_33 * f4t7_33 * f4t8_33 + f2_55 * f4t1_55 * f4t2Dsin_55 * f4t3_55 * f4t4_55 * f4t7_55 * f4t8_55) / rcstackmod);

			// TETA3; t3 = LRACOS (a1 * rhydrodir);
			Ftmp -= stably_normalised(a1 - rcstackdir * cost3) * ((f2_33 * f4t1_33 * f4t2_33 * f4t3Dsin_33 * f4t4_33 * f4t7_33 * f4t8_33 + f2_55 * f4t1_55 * f4t2_55 * f4t3Dsin_55 * f4t4_55 * f4t7_55 * f4t8_55) / rcstackmod);
			Ttmp += stably_normalised(_cross(rcstackdir, a1)) * (-f2_33 * f4t1_33 * f4t2_33 * f4t3Dsin_33 * f4t4_33 * f4t7_33 * f4t8_33 - f2_55 * f4t1_55 * f4t2_55 * f4t3Dsin_55 * f4t4_55 * f4t7_55 * f4t8_55);

			// TETA4; t4 = LRACOS (a3 * b3);
			Ttmp -= stably_normalised(_cross(a3, b3)) * (-f2_33 * f4t1_33 * f4t2_33 * f4t3_33 * f4t4Dsin_33 * f4t7_33 * f4t8_33 - f2_55 * f4t1_55 * f4t2_55 * f4t3_55 * f4t4Dsin_55 * f4t7_55 * f4t8_55);

			// THETA7; t7 = LRACOS (-rcsrackir * b3);
			Ftmp -= stably_normalised(b3 + rcstackdir * cost7) * ((f2_33 * f4t1_33 * f4t2_33 * f4t3_33 * f4t4_33 * f4t7Dsin_33 * f4t8_33 + f2_55 * f4t1_55 * f4t2_55 * f4t3_55 * f4t4_55 * f4t7Dsin_55 * f4t8_55) / rcstackmod);

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			Ftmp -= stably_normalised(a3 - rcstackdir * cost8) * ((f2_33 * f4t1_33 * f4t2_33 * f4t3_33 * f4t4_33 * f4t7_33 * f4t8Dsin_33 + f2_55 * f4t1_55 * f4t2_55 * f4t3_55 * f4t4_55 * f4t7_55 * f4t8Dsin_55) / rcstackmod);
			Ttmp += stably_normalised(_cross(rcstackdir, a3)) * (-f2_33 * f4t1_33 * f4t2_33 * f4t3_33 * f4t4_33 * f4t7_33 * f4t8Dsin_33 - f2_55 * f4t1_55 * f4t2_55 * f4t3_55 * f4t4_55 * f4t7_55 * f4t8Dsin_55);

			Ttmp += _cross(ppos_base, Ftmp);

			Ftmp.w = cstk_energy;
			F += Ftmp;
		}
	}

	// COAXIAL STACKING
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

		c_number f5cosphi3 = 1.f;

		// functions called at their relevant arguments
		c_number f2 = _f2(rstackmod, CXST_F2);
		c_number f4t1 = _f4(t1, CXST_THETA1_T0_OXDNA2, CXST_THETA1_TS, CXST_THETA1_TC, CXST_THETA1_A, CXST_THETA1_B) + _f4_pure_harmonic(t1, CXST_THETA1_SA, CXST_THETA1_SB);
		c_number f4t4 = _f4(t4, CXST_THETA4_T0, CXST_THETA4_TS, CXST_THETA4_TC, CXST_THETA4_A, CXST_THETA4_B);
		c_number f4t5 = _f4(t5, CXST_THETA5_T0, CXST_THETA5_TS, CXST_THETA5_TC, CXST_THETA5_A, CXST_THETA5_B) + _f4(PI - t5, CXST_THETA5_T0, CXST_THETA5_TS, CXST_THETA5_TC, CXST_THETA5_A, CXST_THETA5_B);
		c_number f4t6 = _f4(t6, CXST_THETA6_T0, CXST_THETA6_TS, CXST_THETA6_TC, CXST_THETA6_A, CXST_THETA6_B) + _f4(PI - t6, CXST_THETA6_T0, CXST_THETA6_TS, CXST_THETA6_TC, CXST_THETA6_A, CXST_THETA6_B);

		c_number cxst_energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);

		if(cxst_energy < (c_number) 0) {
			// derivatives called at the relevant arguments
			c_number f2D = _f2D(rstackmod, CXST_F2);
			c_number f4t1D = -_f4D(t1, CXST_THETA1_T0_OXDNA2, CXST_THETA1_TS, CXST_THETA1_TC, CXST_THETA1_A, CXST_THETA1_B) - _f4D_pure_harmonic(t1, CXST_THETA1_SA, CXST_THETA1_SB);
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

			Ttmp += _cross(ppos_stack, Ftmp);
			Ftmp.w = cxst_energy;
			F += Ftmp;
		}
	}
	T += Ttmp;

	// this component stores the energy due to hydrogen bonding
	T.w = tot_Tw;
}

__global__ void DNA3_forces_edge_nonbonded(const c_number4 __restrict__ *poss, const GPU_quat __restrict__ *orientations, c_number4 __restrict__ *forces,
		c_number4 __restrict__ *torques, const edge_bond __restrict__ *edge_list, int n_edges, const int *is_strand_end, bool update_st, CUDAStressTensor *st, const CUDABox *box) {
	if(IND >= n_edges) return;

	c_number4 dF = make_c_number4(0, 0, 0, 0);
	c_number4 dT = make_c_number4(0, 0, 0, 0);

	edge_bond b = edge_list[IND];

	// get info for particle 1
	c_number4 ppos = poss[b.from];
	// particle axes according to Allen's paper
	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[b.from], a1, a2, a3);

	// get info for particle 2
	c_number4 qpos = poss[b.to];
	c_number4 b1, b2, b3;
	get_vectors_from_quat(orientations[b.to], b1, b2, b3);

	c_number4 r = box->minimum_image(ppos, qpos);
	// TODO: DNA3 TO BE UPDATED
	// _DNA3_particle_particle_DNA_interaction(r, ppos, a1, a2, a3, qpos, b1, b2, b3, dF, dT, p_is_end, q_is_end);

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
__global__ void DNA3_forces_edge_bonded(const c_number4 __restrict__ *poss, const GPU_quat __restrict__ *orientations, c_number4 __restrict__ *forces,
		c_number4 __restrict__ *torques, const LR_bonds __restrict__ *bonds, bool use_mbf, c_number mbf_finf, bool update_st, CUDAStressTensor *st) {
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
		uint8_t neigh_n3_type = NO_TYPE; // TODO: DNA3 TO BE UPDATED
		uint8_t neigh_n5_type = NO_TYPE; // TODO: DNA3 TO BE UPDATED
		_DNA3_bonded_part<true>(r, ppos, a1, a2, a3, neigh_n5_type, qpos, b1, b2, b3, neigh_n3_type, dF, T, use_mbf, mbf_finf);
		_update_stress_tensor<true>(p_st, r, dF);
		F += dF;
	}
	if(bs.n5 != P_INVALID) {
		c_number4 qpos = poss[bs.n5];

		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[bs.n5], b1, b2, b3);

		c_number4 r = ppos - qpos;
		c_number4 dF = make_c_number4(0, 0, 0, 0);
		uint8_t neigh_n3_type = NO_TYPE; // TODO: DNA3 TO BE UPDATED
		uint8_t neigh_n5_type = NO_TYPE; // TODO: DNA3 TO BE UPDATED
		_DNA3_bonded_part<false>(r, qpos, b1, b2, b3, neigh_n5_type, ppos, a1, a2, a3, neigh_n3_type, dF, T, use_mbf, mbf_finf);
		_update_stress_tensor<true>(p_st, -r, dF); // -r since r here is defined as r  = ppos - qpos
		F += dF;
	}

	// we can do this because DNA3_forces_edge_nonbonded does not change the reference frame to the torque it calculates
	T = _vectors_transpose_c_number4_product(a1, a2, a3, T);

	forces[IND] = F;
	torques[IND] = T;

	if(update_st) {
		st[IND] += p_st;
	}
}

__global__ void DNA3_forces(const c_number4 __restrict__ *poss, const GPU_quat __restrict__ *orientations, c_number4 __restrict__ *forces,
		c_number4 __restrict__ *torques, const int *matrix_neighs,	const int *number_neighs, const LR_bonds __restrict__ *bonds,
		bool use_mbf, c_number mbf_finf, uint8_t *particle_types, bool update_st, CUDAStressTensor *st, const CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = make_c_number4(0, 0, 0, 0);
	c_number4 ppos = poss[IND];
	LR_bonds pbonds = bonds[IND];
	neigh_types p_neighs(particle_types, pbonds);

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
		uint8_t neigh_n5_type = (pbonds.n5 == P_INVALID) ? NO_TYPE : particle_types[pbonds.n5];
		uint8_t neigh_n3_type = (bonds[pbonds.n3].n3 == P_INVALID) ? NO_TYPE : particle_types[bonds[pbonds.n3].n3];
		_DNA3_bonded_part<true>(r, ppos, a1, a2, a3, neigh_n5_type, qpos, b1, b2, b3, neigh_n3_type, dF, T, use_mbf, mbf_finf);
		_update_stress_tensor<true>(p_st, r, dF);
		F += dF;
	}
	if(pbonds.n5 != P_INVALID) {
		c_number4 qpos = poss[pbonds.n5];
		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[pbonds.n5], b1, b2, b3);

		c_number4 r = ppos - qpos;
		c_number4 dF = make_c_number4(0, 0, 0, 0);
		uint8_t neigh_n5_type = (bonds[pbonds.n5].n5 == P_INVALID) ? NO_TYPE : particle_types[bonds[pbonds.n5].n5];
		uint8_t neigh_n3_type = (pbonds.n3 == P_INVALID) ? NO_TYPE : particle_types[pbonds.n3];
		_DNA3_bonded_part<false>(r, qpos, b1, b2, b3, neigh_n5_type, ppos, a1, a2, a3, neigh_n3_type, dF, T, use_mbf, mbf_finf);
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
			neigh_types q_neighs(particle_types, qbonds);

			c_number4 dF = make_c_number4(0, 0, 0, 0);
			_DNA3_particle_particle_DNA_interaction(r, ppos, a1, a2, a3, qpos, b1, b2, b3, dF, T, p_neighs, q_neighs);
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

__global__ void init_DNA3_strand_ends(int *is_strand_end, const LR_bonds __restrict__ *bonds, int N) {
	if(IND >= N) return;

	LR_bonds pbonds = bonds[IND];
	is_strand_end[IND] = (pbonds.n3 == P_INVALID || pbonds.n5 == P_INVALID);
}

