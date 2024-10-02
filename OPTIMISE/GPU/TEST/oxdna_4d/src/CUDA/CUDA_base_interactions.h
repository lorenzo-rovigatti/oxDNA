
#ifndef CUDA_BASE_INTERACTIONS_H_
#define CUDA_BASE_INTERACTIONS_H_

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
			c_number fmod = 24.f * EXCL_EPS * (lj_part - 2.f*SQR(lj_part)) / rsqr;
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
	int eps_index = 0;
	if(r < MD_F1_RCHIGH[type]) {
		eps_index = 25 * type + n3 * 5 + n5;
		if(r > MD_F1_RHIGH[type]) {
			val = 2.f * MD_F1_BHIGH[type] * (r - MD_F1_RCHIGH[type]);
		}
		else if(r > MD_F1_RLOW[type]) {
			c_number tmp = expf(-(r - MD_F1_R0[type]) * MD_F1_A[type]);
			val = 2.f * (1.f - tmp) * tmp * MD_F1_A[type];
		}
		else if(r > MD_F1_RCLOW[type]) {
			val = 2.f * MD_F1_BLOW[type] * (r - MD_F1_RCLOW[type]);
		}
	}

	return MD_F1_EPS[eps_index] * val;
}


__forceinline__ __device__ c_number _f2(c_number r, int type) {
    c_number val = (c_number) 0.f;
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


__forceinline__ __device__ c_number _f2D(c_number r, int type) {
    c_number val = (c_number) 0.f;
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


__forceinline__ __device__ c_number _f4(c_number t, float t0, float ts, float tc, float a, float b) {
	c_number val = (c_number) 0.f;
	t -= t0;
	if(t < 0) t = -t;

	if(t < tc) {
		if(t > ts) {
			// smoothing
			val = b * SQR(tc - t);
		}
		else val = (c_number) 1.f - a * SQR(t);
	}

	return val;
}


__forceinline__ __device__ c_number _f4Dsin(c_number t, float t0, float ts, float tc, float a, float b) {
	c_number val = (c_number) 0.f;
	c_number tt0 = t - t0;
	// this function is a parabola centered in t0. If tt0 < 0 then the value of the function
	// is the same but the value of its derivative has the opposite sign, so m = -1
	c_number m = copysignf((c_number)1.f, tt0);
	tt0 = copysignf(tt0, (c_number)1.f);

	if(tt0 < tc) {
		c_number sint = sinf(t);
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
#endif
