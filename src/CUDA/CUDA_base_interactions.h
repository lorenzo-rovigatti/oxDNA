
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
#endif
