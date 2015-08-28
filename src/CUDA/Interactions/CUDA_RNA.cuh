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



//__constant__ struct Model RNA_MODEL;



struct CUDAModel {
	float RNA_POS_BACK;
	float RNA_POS_STACK;
	float RNA_POS_BASE;
	// RNA_POS_STACK - RNA_POS_BACK
	float RNA_GAMMA;

	//this is stacking site with your 3-neighbor
	float	RNA_POS_STACK_3_a1;
	float	RNA_POS_STACK_3_a2;
	//this is stacking site for interaction with your 5-neighbor
	float	RNA_POS_STACK_5_a1 ;
	float	RNA_POS_STACK_5_a2 ;
	/*
	 * FENE
	 */
	float RNA_FENE_EPS;
	float RNA_FENE_R0;
	float RNA_FENE_DELTA;
	float RNA_FENE_DELTA2;
	/*
	 * EXCLUDED VOLUME
	 */
	float RNA_EXCL_EPS;
	float RNA_EXCL_S1;
	float RNA_EXCL_S2;
	float RNA_EXCL_S3;
	float RNA_EXCL_S4;
	float RNA_EXCL_R1;
	float RNA_EXCL_R2;
	float RNA_EXCL_R3;
	float RNA_EXCL_R4;
	float RNA_EXCL_B1;
	float RNA_EXCL_B2;
	float RNA_EXCL_B3;
	float RNA_EXCL_B4;
	float RNA_EXCL_RC1;
	float RNA_EXCL_RC2;
	float RNA_EXCL_RC3;
	float RNA_EXCL_RC4;

	/*
	 * HYDROGEN BONDING
	 */
	// radial part
	float RNA_HYDR_EPS;
	float RNA_HYDR_A;
	float RNA_HYDR_RC;
	float RNA_HYDR_R0;
	float RNA_HYDR_BLOW;
	float RNA_HYDR_BHIGH;
	float RNA_HYDR_RLOW;
	float RNA_HYDR_RHIGH;
	float RNA_HYDR_RCLOW;
	float RNA_HYDR_RCHIGH;

	// angular part

	float RNA_HYDR_THETA1_A;
	float RNA_HYDR_THETA1_B;
	float RNA_HYDR_THETA1_T0;
	float RNA_HYDR_THETA1_TS;
	float RNA_HYDR_THETA1_TC;

	float RNA_HYDR_THETA2_A;
	float RNA_HYDR_THETA2_B;
	float RNA_HYDR_THETA2_T0;
	float RNA_HYDR_THETA2_TS;
	float RNA_HYDR_THETA2_TC;

	float RNA_HYDR_THETA3_A;
	float RNA_HYDR_THETA3_B;
	float RNA_HYDR_THETA3_T0;
	float RNA_HYDR_THETA3_TS;
	float RNA_HYDR_THETA3_TC;

	float RNA_HYDR_THETA4_A;
	float RNA_HYDR_THETA4_B;
	float RNA_HYDR_THETA4_T0;
	float RNA_HYDR_THETA4_TS;
	float RNA_HYDR_THETA4_TC;

	float RNA_HYDR_THETA7_A;
	float RNA_HYDR_THETA7_B;
	float RNA_HYDR_THETA7_T0;
	float RNA_HYDR_THETA7_TS;
	float RNA_HYDR_THETA7_TC;

	float RNA_HYDR_THETA8_A;
	float RNA_HYDR_THETA8_B;
	float RNA_HYDR_THETA8_T0;
	float RNA_HYDR_THETA8_TS;
	float RNA_HYDR_THETA8_TC;

	/*
	 * STACKING
	 */
	// radial part
	float RNA_STCK_BASE_EPS;
	float RNA_STCK_FACT_EPS;
	float RNA_STCK_A;
	float RNA_STCK_RC;
	float RNA_STCK_R0;
	float RNA_STCK_BLOW;
	float RNA_STCK_BHIGH;
	float RNA_STCK_RLOW;
	float RNA_STCK_RHIGH;
	float RNA_STCK_RCLOW;
	float RNA_STCK_RCHIGH;

	// angular part, theta
	float RNA_STCK_THETA4_A;
	float RNA_STCK_THETA4_B;
	float RNA_STCK_THETA4_T0;
	float RNA_STCK_THETA4_TS;
	float RNA_STCK_THETA4_TC;
	float RNA_STCK_THETA5_A;
	float RNA_STCK_THETA5_B;
	float RNA_STCK_THETA5_T0;
	float RNA_STCK_THETA5_TS;
	float RNA_STCK_THETA5_TC;
	float RNA_STCK_THETA6_A;
	float RNA_STCK_THETA6_B;
	float RNA_STCK_THETA6_T0;
	float RNA_STCK_THETA6_TS;
	float RNA_STCK_THETA6_TC;

	float STCK_THETAB1_A;
	float STCK_THETAB1_B;
	float STCK_THETAB1_T0;
	float STCK_THETAB1_TS;
	float STCK_THETAB1_TC;

	float STCK_THETAB2_A;
	float STCK_THETAB2_B;
	float STCK_THETAB2_T0;
	float STCK_THETAB2_TS;
	float STCK_THETAB2_TC;




	// angular part, phi
	float RNA_STCK_PHI1_A;
	float RNA_STCK_PHI1_B;
	float RNA_STCK_PHI1_XC;
	float RNA_STCK_PHI1_XS;
	float RNA_STCK_PHI2_A;
	float RNA_STCK_PHI2_B;
	float RNA_STCK_PHI2_XC;
	float RNA_STCK_PHI2_XS;

	/*
	 *  CROSS STACKING
	 */
	// radial part
	float RNA_CRST_R0;
	float RNA_CRST_RC;
	float RNA_CRST_K;
	float RNA_CRST_BLOW;
	float RNA_CRST_RLOW;
	float RNA_CRST_RCLOW;
	float RNA_CRST_BHIGH;
	float RNA_CRST_RHIGH;
	float RNA_CRST_RCHIGH;
	// angular part;
	float RNA_CRST_THETA1_A;
	float RNA_CRST_THETA1_B;
	float RNA_CRST_THETA1_T0;
	float RNA_CRST_THETA1_TS;
	float RNA_CRST_THETA1_TC;
	float RNA_CRST_THETA2_A;
	float RNA_CRST_THETA2_B;
	float RNA_CRST_THETA2_T0;
	float RNA_CRST_THETA2_TS;
	float RNA_CRST_THETA2_TC;
	float RNA_CRST_THETA3_A;
	float RNA_CRST_THETA3_B;
	float RNA_CRST_THETA3_T0;
	float RNA_CRST_THETA3_TS;
	float RNA_CRST_THETA3_TC;
	float RNA_CRST_THETA4_A;
	float RNA_CRST_THETA4_B;
	float RNA_CRST_THETA4_T0;
	float RNA_CRST_THETA4_TS;
	float RNA_CRST_THETA4_TC;
	float RNA_CRST_THETA7_A;
	float RNA_CRST_THETA7_B;
	float RNA_CRST_THETA7_T0;
	float RNA_CRST_THETA7_TS;
	float RNA_CRST_THETA7_TC;
	float RNA_CRST_THETA8_A;
	float RNA_CRST_THETA8_B;
	float RNA_CRST_THETA8_T0;
	float RNA_CRST_THETA8_TS;
	float RNA_CRST_THETA8_TC;

	/*
	 *  COAXIAL STACKING
	 */
	// radial part
	float RNA_CXST_R0;
	float RNA_CXST_RC;
	float RNA_CXST_K;
	float RNA_CXST_BLOW;
	float RNA_CXST_RLOW;
	float RNA_CXST_RCLOW;
	float RNA_CXST_BHIGH;
	float RNA_CXST_RHIGH;
	float RNA_CXST_RCHIGH;
	// angular part;

	float RNA_CXST_THETA1_A;
	float RNA_CXST_THETA1_B;
	float RNA_CXST_THETA1_T0;
	float RNA_CXST_THETA1_TS;
	float RNA_CXST_THETA1_TC;
	float RNA_CXST_THETA4_A;
	float RNA_CXST_THETA4_B;
	float RNA_CXST_THETA4_T0;
	float RNA_CXST_THETA4_TS;
	float RNA_CXST_THETA4_TC;
	float RNA_CXST_THETA5_A;
	float RNA_CXST_THETA5_B;
	float RNA_CXST_THETA5_T0;
	float RNA_CXST_THETA5_TS;
	float RNA_CXST_THETA5_TC;
	float RNA_CXST_THETA6_A;
	float RNA_CXST_THETA6_B;
	float RNA_CXST_THETA6_T0;
	float RNA_CXST_THETA6_TS;
	float RNA_CXST_THETA6_TC;
	// angular part, phi3 and phi4

	float RNA_CXST_PHI3_A;
	float RNA_CXST_PHI3_B;
	float RNA_CXST_PHI3_XC;
	float RNA_CXST_PHI3_XS;
	float RNA_CXST_PHI4_A;
	float RNA_CXST_PHI4_B;
	float RNA_CXST_PHI4_XC;
	float RNA_CXST_PHI4_XS;

//new stuff
	//float INCLINATION;
    //float BASE_SIZE_R;
    //float ROT_SUGAR;
    float p3_x, p3_y, p3_z, p5_x, p5_y, p5_z;
    float RNA_POS_BACK_a1, RNA_POS_BACK_a2, RNA_POS_BACK_a3;
};


__constant__ CUDAModel rnamodel;





#include "../cuda_utils/CUDA_lr_common.cuh"

template<typename number, typename number4>
__forceinline__ __device__ void _excluded_volume(const number4 &r, number4 &F, number sigma, number rstar, number b, number rc) {
	number rsqr = CUDA_DOT(r, r);

	F.x = F.y = F.z = F.w = (number) 0.f;
	if(rsqr < SQR(rc)) {
		if(rsqr > SQR(rstar)) {
			number rmod = sqrt(rsqr);
			number rrc = rmod - rc;
			number fmod = 2.f * rnamodel.RNA_EXCL_EPS * b * rrc / rmod;
			F.x = r.x * fmod;
			F.y = r.y * fmod;
			F.z = r.z * fmod;
			F.w = rnamodel.RNA_EXCL_EPS * b * SQR(rrc);
		}
		else {
			number lj_part = CUB(SQR(sigma)/rsqr);
			number fmod = 24.f * rnamodel.RNA_EXCL_EPS * (lj_part - 2.f*SQR(lj_part)) / rsqr;
			F.x = r.x * fmod;
			F.y = r.y * fmod;
			F.z = r.z * fmod;
			F.w = 4.f * rnamodel.RNA_EXCL_EPS * (SQR(lj_part) - lj_part);
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
__forceinline__ __device__ number _fX(number r, int type, int n3, int n5) {
	number val = (number) 0.f;
	if(r < MD_F1_RCHIGH[type]) {
		int eps_index = 25 * type + n3 * 5 + n5;
		if(r > MD_F1_RHIGH[type]) {
			val = MD_F1_EPS[eps_index] * MD_F1_BHIGH[type] * SQR(r - MD_F1_RCHIGH[type]);
		}
		else if(r > MD_F1_R0[type]) {
			number tmp = 1.f - expf(-(r - MD_F1_R0[type]) * MD_F1_A[type]);
			val = MD_F1_EPS[eps_index] * SQR(tmp) - MD_F1_SHIFT[eps_index];
		}
		else  {
		     val = - MD_F1_SHIFT[eps_index];
			//val = MD_F1_EPS[eps_index] * MD_F1_BLOW[type] * SQR(r - MD_F1_RCLOW[type]);
		}
	}

	return val;
}

template<typename number>
__forceinline__ __device__ number _fXD(number r, int type, int n3, int n5) {
	number val = (number) 0.f;
	int eps_index = 0;
	if(r < MD_F1_RCHIGH[type]) {
		eps_index = 25 * type + n3 * 5 + n5;
		if(r > MD_F1_RHIGH[type]) {
			val = 2.f * MD_F1_BHIGH[type] * (r - MD_F1_RCHIGH[type]);
		}
		else if(r > MD_F1_R0[type]) {
			number tmp = expf(-(r - MD_F1_R0[type]) * MD_F1_A[type]);
			val = 2.f * (1.f - tmp) * tmp * MD_F1_A[type];
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
/*
	if(type == 1)
		printf("DEVICE _f2D: %f,  Running with %f %f %f %f %f %f %f\n",val,MD_F2_RCHIGH[type],MD_F2_RHIGH[type],MD_F2_BHIGH[type],MD_F2_BHIGH[type],MD_F2_RCLOW[type],MD_F2_BLOW[type], MD_F2_K[type]);
*/
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
	_excluded_volume(rcenter, Ftmp, rnamodel.RNA_EXCL_S2, rnamodel.RNA_EXCL_R2, rnamodel.RNA_EXCL_B2, rnamodel.RNA_EXCL_RC2);
	number4 torquep1 = (qIsN3) ? _cross<number, number4>(n5pos_base, Ftmp) : _cross<number, number4>(n3pos_base, Ftmp);
	F += Ftmp;

	// n5-BASE vs. n3-BACK
	rcenter = r + n3pos_back - n5pos_base;
	_excluded_volume(rcenter, Ftmp, rnamodel.RNA_EXCL_S3, rnamodel.RNA_EXCL_R3, rnamodel.RNA_EXCL_B3, rnamodel.RNA_EXCL_RC3);
	number4 torquep2 = (qIsN3) ? _cross<number, number4>(n5pos_base, Ftmp) : _cross<number, number4>(n3pos_back, Ftmp);
	F += Ftmp;

	// n5-BACK vs. n3-BASE
	rcenter = r + n3pos_base - n5pos_back;
	_excluded_volume(rcenter, Ftmp, rnamodel.RNA_EXCL_S4, rnamodel.RNA_EXCL_R4, rnamodel.RNA_EXCL_B4, rnamodel.RNA_EXCL_RC4);
	number4 torquep3 = (qIsN3) ? _cross<number, number4>(n5pos_back, Ftmp) : _cross<number, number4>(n3pos_base, Ftmp);
	F += Ftmp;

	T += torquep1 + torquep2 + torquep3;
}

template <typename number, typename number4, bool qIsN3> __device__
void _bonded_part(number4 &n5pos, number4 &n5x, number4 &n5y, number4 &n5z,number4 &n3pos,
			     number4 &n3x, number4 &n3y, number4 &n3z, number4 &F, number4 &T, bool average) {

	//printf("Hello from bonded part function \n");
	int n3type = get_particle_type<number, number4>(n3pos);
	int n5type = get_particle_type<number, number4>(n5pos);

	number4 r = make_number4<number, number4>(n3pos.x - n5pos.x, n3pos.y - n5pos.y, n3pos.z - n5pos.z, (number) 0);

	number4 n5pos_back;
	number4 n3pos_back;
	//if(grooving) n5pos_back = n5x * POS_MM_BACK1 + n5y * POS_MM_BACK2;
	//else n5pos_back = n5x * rnamodel.RNA_POS_BACK;

	n3pos_back = n3x * rnamodel.RNA_POS_BACK_a1 + n3y * rnamodel.RNA_POS_BACK_a2 + n3z * rnamodel.RNA_POS_BACK_a3;
	n5pos_back = n5x * rnamodel.RNA_POS_BACK_a1 + n5y * rnamodel.RNA_POS_BACK_a2 + n5z * rnamodel.RNA_POS_BACK_a3;


	number4 n5pos_base = n5x * rnamodel.RNA_POS_BASE;




	//if(grooving) n3pos_back = n3x * POS_MM_BACK1 + n3y * POS_MM_BACK2;
	//else n3pos_back = n3x * rnamodel.RNA_POS_BACK;
	number4 n3pos_base = n3x * rnamodel.RNA_POS_BASE;
	number4 n3pos_stack_5 = n3x * rnamodel.RNA_POS_STACK_5_a1   + n3y * rnamodel.RNA_POS_STACK_5_a2;  //n3x * rnamodel.RNA_POS_STACK;
	//number4 n3pos_stack_3 = n3x * rnamodel.RNA_POS_STACK_3_a1   + n3y * rnamodel.RNA_POS_STACK_3_a2;  //n3x * rnamodel.RNA_POS_STACK;

	//number4 n5pos_stack_5 = n5x * rnamodel.RNA_POS_STACK_5_a1   + n5y * rnamodel.RNA_POS_STACK_5_a2;  //n3x * rnamodel.RNA_POS_STACK;
	number4 n5pos_stack_3 = n5x * rnamodel.RNA_POS_STACK_3_a1   + n5y * rnamodel.RNA_POS_STACK_3_a2;  //n3x * rnamodel.RNA_POS_STACK;


	number4 rback = r + n3pos_back - n5pos_back;

	number rbackmod = _module<number, number4>(rback);

	number4 rbackdir = make_number4<number,number4>(rback.x / rbackmod, rback.y / rbackmod, rback.z / rbackmod, 0);

	number rbackr0 = rbackmod - rnamodel.RNA_FENE_R0;

	number4 Ftmp = rback * ((rnamodel.RNA_FENE_EPS * rbackr0  / (rnamodel.RNA_FENE_DELTA2 - SQR(rbackr0))) / rbackmod);
	Ftmp.w = -rnamodel.RNA_FENE_EPS * ((number)0.5f) * logf(1 - SQR(rbackr0) / rnamodel.RNA_FENE_DELTA2);

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
	number4 rstack =  r + n3pos_stack_5 - n5pos_stack_3;

	number rstackmod = _module<number, number4>(rstack);
	number4 rstackdir = make_number4<number, number4>(rstack.x / rstackmod, rstack.y / rstackmod, rstack.z / rstackmod, 0);
	// This is the position the backbone would have with major-minor grooves the same width.
	// We need to do this to implement different major-minor groove widths because rback is
	// used as a reference point for things that have nothing to do with the actual backbone
	// position (in this case, the stacking interaction).
	//number4 rbackref = r + n3x * rnamodel.RNA_POS_BACK - n5x * rnamodel.RNA_POS_BACK;
	//number rbackrefmod = _module<number, number4>(rbackref);

	number t4 = CUDA_LRACOS(CUDA_DOT(n3z, n5z));

	number t5 = CUDA_LRACOS(CUDA_DOT(n5z, rstackdir));
	number t6 = CUDA_LRACOS(-CUDA_DOT(n3z, rstackdir));


	number cosphi1 = CUDA_DOT(n5y, rbackdir); // / rbackrefmod;
	number cosphi2 = CUDA_DOT(n3y, rbackdir); // / rbackrefmod;


	//number costB1 = -rback_back_dir * p->int_centers[RNANucleotide<number>::BBVECTOR_3];
	//number costB2 = -rback_back_dir * q->int_centers[RNANucleotide<number>::BBVECTOR_5];

	number4 n3bbvector_5 = (n3x * rnamodel.p5_x + n3y * rnamodel.p5_y +  n3z * rnamodel.p5_z);
	number4 n5bbvector_3 = (n5x * rnamodel.p3_x + n5y * rnamodel.p3_y +  n5z * rnamodel.p3_z);


	number costB1 = -CUDA_DOT( rbackdir,n5bbvector_3 ) ;
	number costB2 = -CUDA_DOT( rbackdir,n3bbvector_5 )   ;
	number tB1 = CUDA_LRACOS( costB1 );
	number tB2 = CUDA_LRACOS( costB2 );


	// functions
	number f1 = _f1(rstackmod, STCK_F1, n3type, n5type);
	//number f4t4 = _f4(t4, rnamodel.RNA_STCK_THETA4_T0, rnamodel.RNA_STCK_THETA4_TS, rnamodel.RNA_STCK_THETA4_TC, rnamodel.RNA_STCK_THETA4_A, rnamodel.RNA_STCK_THETA4_B);
	number f4t5 = _f4(PI - t5, rnamodel.RNA_STCK_THETA5_T0, rnamodel.RNA_STCK_THETA5_TS, rnamodel.RNA_STCK_THETA5_TC, rnamodel.RNA_STCK_THETA5_A, rnamodel.RNA_STCK_THETA5_B);
	number f4t6 = _f4(t6, rnamodel.RNA_STCK_THETA6_T0, rnamodel.RNA_STCK_THETA6_TS, rnamodel.RNA_STCK_THETA6_TC, rnamodel.RNA_STCK_THETA6_A, rnamodel.RNA_STCK_THETA6_B);

	number f4tB1 = _f4(tB1, rnamodel.STCK_THETAB1_T0, rnamodel.STCK_THETAB1_TS, rnamodel.STCK_THETAB1_TC, rnamodel.STCK_THETAB1_A, rnamodel.STCK_THETAB1_B);
	number f4tB2 = _f4(tB2, rnamodel.STCK_THETAB2_T0, rnamodel.STCK_THETAB2_TS, rnamodel.STCK_THETAB2_TC, rnamodel.STCK_THETAB2_A, rnamodel.STCK_THETAB2_B);



	number f5phi1 = _f5(cosphi1, STCK_F5_PHI1);
	number f5phi2 = _f5(cosphi2, STCK_F5_PHI2);

	number energy = f1 * f4tB1 * f4tB2 * f4t5 * f4t6 * f5phi1 * f5phi2;


	//number4 pre_Ttmp = Ttmp;
	//number4 pre_Ftmp = Ftmp;

	if(energy != (number) 0) {
		// and their derivatives
		number f1D = _f1D(rstackmod, STCK_F1, n3type, n5type);
		//number f4t4Dsin = _f4Dsin(t4, rnamodel.RNA_STCK_THETA4_T0, rnamodel.RNA_STCK_THETA4_TS, rnamodel.RNA_STCK_THETA4_TC, rnamodel.RNA_STCK_THETA4_A, rnamodel.RNA_STCK_THETA4_B);
		number f4t5Dsin = _f4Dsin(PI - t5, rnamodel.RNA_STCK_THETA5_T0, rnamodel.RNA_STCK_THETA5_TS, rnamodel.RNA_STCK_THETA5_TC, rnamodel.RNA_STCK_THETA5_A, rnamodel.RNA_STCK_THETA5_B);
		number f4t6Dsin = _f4Dsin(t6, rnamodel.RNA_STCK_THETA6_T0, rnamodel.RNA_STCK_THETA6_TS, rnamodel.RNA_STCK_THETA6_TC, rnamodel.RNA_STCK_THETA6_A, rnamodel.RNA_STCK_THETA6_B);

		number f5phi1D = _f5D(cosphi1, STCK_F5_PHI1);
		number f5phi2D = _f5D(cosphi2, STCK_F5_PHI2);

		number f4tB1Dsin = _f4Dsin(tB1, rnamodel.STCK_THETAB1_T0, rnamodel.STCK_THETAB1_TS, rnamodel.STCK_THETAB1_TC, rnamodel.STCK_THETAB1_A, rnamodel.STCK_THETAB1_B);
		number f4tB2Dsin = _f4Dsin(tB2, rnamodel.STCK_THETAB2_T0, rnamodel.STCK_THETAB2_TS, rnamodel.STCK_THETAB2_TC, rnamodel.STCK_THETAB2_A, rnamodel.STCK_THETAB2_B);



		// RADIAL
		Ftmp = -rstackdir * (energy * f1D / f1);

		//if(grooving) printf("Amydevice: (qisN3 = %d) calculated Ftmp=(%f,%f,%f) \n",qIsN3,Ftmp.x,Ftmp.y,Ftmp.z);
		// THETA 5
		Ftmp -= (n5z - cosf(t5) * rstackdir) * (f1 * f4t5Dsin * f4t6 * f5phi1 * f5phi2 * f4tB1 * f4tB2 / rstackmod); //(n5z - cosf(t5) * rstackdir) * (energy * f4t5Dsin / (f4t5 * rstackmod));
	  // THETA 6
		Ftmp -= (n3z + cosf(t6) * rstackdir) * (energy * f4t6Dsin / (f4t6 * rstackmod));
		//if(grooving) printf("mydevice: (qisN3 = %d) calculated Ftmp=(%f,%f,%f) \n ",qIsN3,Ftmp.x,Ftmp.y,Ftmp.z);

		number4 Fposback  = -( n5bbvector_3 + rbackdir * costB1) * (f1 * f4t5 * f4t6 * f5phi1 * f5phi2 * f4tB1Dsin * f4tB2 / rbackmod );
		Fposback -=  (n3bbvector_5 + rbackdir * costB2) * (f1 * f4t5 * f4t6 * f5phi1 * f5phi2 * f4tB1 * f4tB2Dsin / rbackmod );


  //force acting on the backbone site

		// COS PHI 1
		// here particle p is referred to using the a while particle q is referred with the b
		/*
		number ra2 = CUDA_DOT(rstackdir, n5y);
		number ra1 = CUDA_DOT(rstackdir, n5x);
		number rb1 = CUDA_DOT(rstackdir, n3x);
		number a2b1 = CUDA_DOT(n5y, n3x);
		number dcosphi1dr = (SQR(rstackmod)*ra2 - ra2*SQR(rbackrefmod) - rstackmod*(a2b1 + ra2*(-ra1 + rb1))*rnamodel.RNA_GAMMA + a2b1*(-ra1 + rb1)*SQR(rnamodel.RNA_GAMMA))/(SQR(rbackrefmod)*rbackrefmod);
		number dcosphi1dra1 = rstackmod*rnamodel.RNA_GAMMA*(rstackmod*ra2 - a2b1*rnamodel.RNA_GAMMA)/(SQR(rbackrefmod)*rbackrefmod);
		number dcosphi1dra2 = -rstackmod / rbackrefmod;
		number dcosphi1drb1 = -(rstackmod*rnamodel.RNA_GAMMA*(rstackmod*ra2 - a2b1*rnamodel.RNA_GAMMA))/(SQR(rbackrefmod)*rbackrefmod);
		number dcosphi1da1b1 = SQR(rnamodel.RNA_GAMMA)*(-rstackmod*ra2 + a2b1*rnamodel.RNA_GAMMA)/(SQR(rbackrefmod)*rbackrefmod);
		number dcosphi1da2b1 = rnamodel.RNA_GAMMA / rbackrefmod;
        */
		//number force_part_phi1 = energy * f5phi1D / f5phi1;

		number force_part_phi1 = -f1 * f4t5 * f4t6 * f5phi1D * f5phi2 * f4tB1 * f4tB2;
		Fposback += (n5y - rback * (cosphi1 / rbackmod) ) * (force_part_phi1 / rbackmod) ;


		/*
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
		number dcosphi2dr = ((rstackmod*ra2 + a2b1*rnamodel.RNA_GAMMA)*(rstackmod + (rb1 - ra1)*rnamodel.RNA_GAMMA) - ra2*SQR(rbackrefmod))/(SQR(rbackrefmod)*rbackrefmod);
		number dcosphi2dra1 = -rstackmod*rnamodel.RNA_GAMMA*(rstackmod*ra2 + a2b1*rnamodel.RNA_GAMMA)/(SQR(rbackrefmod)*rbackrefmod);
		number dcosphi2dra2 = -rstackmod / rbackrefmod;
		number dcosphi2drb1 = (rstackmod*rnamodel.RNA_GAMMA*(rstackmod*ra2 + a2b1*rnamodel.RNA_GAMMA))/(SQR(rbackrefmod)*rbackrefmod);
		number dcosphi2da1b1 = -SQR(rnamodel.RNA_GAMMA)*(rstackmod*ra2 + a2b1*rnamodel.RNA_GAMMA)/(SQR(rbackrefmod)*rbackrefmod);
		number dcosphi2da2b1 = -rnamodel.RNA_GAMMA / rbackrefmod;

		number force_part_phi2 = energy * f5phi2D / f5phi2;

		Ftmp -= (rstackdir * dcosphi2dr +
			    ((n3y - rstackdir * ra2) * dcosphi2dra2 +
				(n3x - rstackdir * ra1) * dcosphi2dra1 +
				(n5x - rstackdir * rb1) * dcosphi2drb1) / rstackmod) * force_part_phi2;

        */

		number force_part_phi2 = -f1  * f4t5 * f4t6 * f5phi1 * f5phi2D;
		Fposback += (n3y - rback * (cosphi2 / rbackmod) ) * (force_part_phi2 / rbackmod) ;
		//if(grooving) printf("mydevice: (qisN3 = %d) calculated Fposback=(%f,%f,%f) \n ",qIsN3,Fposback.x,Fposback.y,Fposback.z);


		// TORQUE
		if(qIsN3) {
			Ttmp = -_cross<number, number4>(n5pos_stack_3, Ftmp);
			Ttmp -= _cross<number,number4>(n5pos_back, Fposback);



			// THETA 5
			Ttmp += _cross<number,number4>(rstackdir,n5z) * (f1  * f4t5Dsin * f4t6 * f5phi1 * f5phi2 * f4tB1 * f4tB2);

		    //thetaB1
			Ttmp +=  _cross<number,number4>(rbackdir,n5bbvector_3) * (f1  * f4t5 * f4t6 * f5phi1 * f5phi2 * f4tB1Dsin * f4tB2);

			Ttmp += _cross<number,number4>(n5y,rbackdir) * force_part_phi1;



		}
		else {

			Ttmp = _cross<number, number4>(n3pos_stack_5, Ftmp);
			Ttmp += _cross<number,number4>(n3pos_back, Fposback);
			Ttmp += _cross<number,number4>(n3y,rbackdir) * force_part_phi2;

			// THETA 6
			//printf("Entered function with qIsN3 = %d\n",(int)qIsN3);
			Ttmp += _cross<number,number4>(rstackdir,n3z) *f1 * f4t5 * f4t6Dsin * f5phi1 * f5phi2 * f4tB1 * f4tB2;

			//theta B2
			Ttmp += _cross<number,number4>(rbackdir,n3bbvector_5) * ( f1 * f4t5 * f4t6 * f5phi1 * f5phi2 * f4tB1 * f4tB2Dsin );

		}


		Ftmp.w = energy;


		//printf("Approaching decision loop with qIsN3 = %d\n",(int)qIsN3);
		if(qIsN3) {

            Ftmp += Fposback;


		    T += Ttmp;
			F -= Ftmp;

		}
		else {
			Ftmp += Fposback;



			T += Ttmp;
			F += Ftmp ;
		}


	}

   //Ftmp = F;
	//if(grooving)
	//  printf("This is bonded function that calculated F=(%f,%f,%f) and T=(%f,%f,%f)  and qIsN3 is %d \n",Ftmp.x,Ftmp.y,Ftmp.z,Ttmp.x,Ttmp.y,Ttmp.z,(int)qIsN3);


	//printf("Goodbye from bonded part function \n");
}

template <typename number, typename number4> __device__ void _particle_particle_interaction(number4 ppos, number4 a1, number4 a2, number4 a3, number4 qpos, number4 b1, number4 b2, number4 b3, number4 &F, number4 &T, bool average, bool use_debye_huckel, bool mismatch_repulsion, LR_bonds pbonds, LR_bonds qbonds) {
	int ptype = get_particle_type<number, number4>(ppos);
	int qtype = get_particle_type<number, number4>(qpos);
	int pbtype = get_particle_btype<number, number4>(ppos);
	int qbtype = get_particle_btype<number, number4>(qpos);
	int int_type = pbtype + qbtype;

	number4 r = minimum_image<number, number4>(ppos, qpos);

	number4 ppos_back = a1 * rnamodel.RNA_POS_BACK_a1 + a2 * rnamodel.RNA_POS_BACK_a2 + a3 * rnamodel.RNA_POS_BACK_a3;
	//if(grooving) ppos_back = POS_MM_BACK1 * a1 + POS_MM_BACK2 * a2;
  	//else ppos_back = rnamodel.RNA_POS_BACK * a1;

	number4 ppos_base = rnamodel.RNA_POS_BASE * a1;
	number4 ppos_stack = rnamodel.RNA_POS_STACK * a1;

	number4 qpos_back = b1 * rnamodel.RNA_POS_BACK_a1 + b2 * rnamodel.RNA_POS_BACK_a2 + b3 * rnamodel.RNA_POS_BACK_a3;
	//if(grooving) qpos_back = POS_MM_BACK1 * b1 + POS_MM_BACK2 * b2;
    //else qpos_back = rnamodel.RNA_POS_BACK * b1;

	number4 qpos_base = rnamodel.RNA_POS_BASE * b1;
	number4 qpos_stack = rnamodel.RNA_POS_STACK * b1;

	number old_Tw = T.w;

	// excluded volume
	// BACK-BACK
	number4 Ftmp = make_number4<number, number4>(0, 0, 0, 0);
	number4 rbackbone = r + qpos_back - ppos_back;
	_excluded_volume(rbackbone, Ftmp, rnamodel.RNA_EXCL_S1, rnamodel.RNA_EXCL_R1, rnamodel.RNA_EXCL_B1, rnamodel.RNA_EXCL_RC1);
	number4 Ttmp = _cross<number, number4>(ppos_back, Ftmp);
	_bonded_excluded_volume<number, number4, true>(r, qpos_base, qpos_back, ppos_base, ppos_back, Ftmp, Ttmp);

	F += Ftmp;

	// HYDROGEN BONDING
	int is_pair = int(int_type == 3);
	if(!average) //allow for wobble bp
	{
			if( int_type == 4 && ( (qtype == N_T && ptype == N_G ) || (qtype == N_G && ptype == N_T)  ))
				is_pair = 1;
	}

	number hb_energy = (number) 0;
	number4 rhydro = r + qpos_base - ppos_base;
	number rhydromodsqr = CUDA_DOT(rhydro, rhydro);
	if( (is_pair || mismatch_repulsion) && SQR(rnamodel.RNA_HYDR_RCLOW) < rhydromodsqr && rhydromodsqr < SQR(rnamodel.RNA_HYDR_RCHIGH)) {
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
		number f1 = hb_multi * _f1(rhydromod, RNA_HYDR_F1, ptype, qtype);
		if(mismatch_repulsion && !is_pair)
		{
		    f1 =  _fX(rhydromod, RNA_HYDR_F1, 0, 0);
		}
		number f4t1 = _f4(t1, rnamodel.RNA_HYDR_THETA1_T0, rnamodel.RNA_HYDR_THETA1_TS, rnamodel.RNA_HYDR_THETA1_TC, rnamodel.RNA_HYDR_THETA1_A, rnamodel.RNA_HYDR_THETA1_B);
		number f4t2 = _f4(t2, rnamodel.RNA_HYDR_THETA2_T0, rnamodel.RNA_HYDR_THETA2_TS, rnamodel.RNA_HYDR_THETA2_TC, rnamodel.RNA_HYDR_THETA2_A, rnamodel.RNA_HYDR_THETA2_B);
		number f4t3 = _f4(t3, rnamodel.RNA_HYDR_THETA3_T0, rnamodel.RNA_HYDR_THETA3_TS, rnamodel.RNA_HYDR_THETA3_TC, rnamodel.RNA_HYDR_THETA3_A, rnamodel.RNA_HYDR_THETA3_B);
		number f4t4 = _f4(t4, rnamodel.RNA_HYDR_THETA4_T0, rnamodel.RNA_HYDR_THETA4_TS, rnamodel.RNA_HYDR_THETA4_TC, rnamodel.RNA_HYDR_THETA4_A, rnamodel.RNA_HYDR_THETA4_B);
		number f4t7 = _f4(t7, rnamodel.RNA_HYDR_THETA7_T0, rnamodel.RNA_HYDR_THETA7_TS, rnamodel.RNA_HYDR_THETA7_TC, rnamodel.RNA_HYDR_THETA7_A, rnamodel.RNA_HYDR_THETA7_B);
		number f4t8 = _f4(t8, rnamodel.RNA_HYDR_THETA8_T0, rnamodel.RNA_HYDR_THETA8_TS, rnamodel.RNA_HYDR_THETA8_TC, rnamodel.RNA_HYDR_THETA8_A, rnamodel.RNA_HYDR_THETA8_B);

		hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

		if(hb_energy < (number) 0  ||  (hb_energy > 0 && mismatch_repulsion && !is_pair)) {
			// derivatives called at the relevant arguments
			number f1D = hb_multi * _f1D(rhydromod, HYDR_F1, ptype, qtype);
			if(mismatch_repulsion && !is_pair)
			{
			   f1D =  _fXD(rhydromod, RNA_HYDR_F1, 0, 0);
			}
			number f4t1Dsin = - _f4Dsin(t1, rnamodel.RNA_HYDR_THETA1_T0, rnamodel.RNA_HYDR_THETA1_TS, rnamodel.RNA_HYDR_THETA1_TC, rnamodel.RNA_HYDR_THETA1_A, rnamodel.RNA_HYDR_THETA1_B);
			number f4t2Dsin = - _f4Dsin(t2, rnamodel.RNA_HYDR_THETA2_T0, rnamodel.RNA_HYDR_THETA2_TS, rnamodel.RNA_HYDR_THETA2_TC, rnamodel.RNA_HYDR_THETA2_A, rnamodel.RNA_HYDR_THETA2_B);
			number f4t3Dsin = _f4Dsin(t3, rnamodel.RNA_HYDR_THETA3_T0, rnamodel.RNA_HYDR_THETA3_TS, rnamodel.RNA_HYDR_THETA3_TC, rnamodel.RNA_HYDR_THETA3_A, rnamodel.RNA_HYDR_THETA3_B);
			number f4t4Dsin = _f4Dsin(t4, rnamodel.RNA_HYDR_THETA4_T0, rnamodel.RNA_HYDR_THETA4_TS, rnamodel.RNA_HYDR_THETA4_TC, rnamodel.RNA_HYDR_THETA4_A, rnamodel.RNA_HYDR_THETA4_B);
			number f4t7Dsin = - _f4Dsin(t7, rnamodel.RNA_HYDR_THETA7_T0, rnamodel.RNA_HYDR_THETA7_TS, rnamodel.RNA_HYDR_THETA7_TC, rnamodel.RNA_HYDR_THETA7_A, rnamodel.RNA_HYDR_THETA7_B);
			number f4t8Dsin = _f4Dsin(t8, rnamodel.RNA_HYDR_THETA8_T0, rnamodel.RNA_HYDR_THETA8_TS, rnamodel.RNA_HYDR_THETA8_TC, rnamodel.RNA_HYDR_THETA8_A, rnamodel.RNA_HYDR_THETA8_B);

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
	if(SQR(rnamodel.RNA_CRST_RCLOW) < rcstackmodsqr && rcstackmodsqr < SQR(rnamodel.RNA_CRST_RCHIGH)) {
	  	number rcstackmod = sqrtf(rcstackmodsqr);
	  	number4 rcstackdir = rcstack / rcstackmod;

		// angles involved in the CSTCK interaction
		number t1 = CUDA_LRACOS (-CUDA_DOT(a1, b1));

		number t2 = CUDA_LRACOS (-CUDA_DOT(b1, rcstackdir));
		number t3 = CUDA_LRACOS ( CUDA_DOT(a1, rcstackdir));

		//number t4 = CUDA_LRACOS ( CUDA_DOT(a3, b3));

		number t7 = CUDA_LRACOS (-CUDA_DOT(rcstackdir, b3));
		number t8 = CUDA_LRACOS ( CUDA_DOT(rcstackdir, a3));

	 	 // functions called at their relevant arguments
		number f2 = _f2(rcstackmod, RNA_CRST_F2);
		number f4t1 = _f4(t1, rnamodel.RNA_CRST_THETA1_T0, rnamodel.RNA_CRST_THETA1_TS, rnamodel.RNA_CRST_THETA1_TC, rnamodel.RNA_CRST_THETA1_A, rnamodel.RNA_CRST_THETA1_B);
		number f4t2 = _f4(t2, rnamodel.RNA_CRST_THETA2_T0, rnamodel.RNA_CRST_THETA2_TS, rnamodel.RNA_CRST_THETA2_TC, rnamodel.RNA_CRST_THETA2_A, rnamodel.RNA_CRST_THETA2_B);
		number f4t3 = _f4(t3, rnamodel.RNA_CRST_THETA3_T0, rnamodel.RNA_CRST_THETA3_TS, rnamodel.RNA_CRST_THETA3_TC, rnamodel.RNA_CRST_THETA3_A, rnamodel.RNA_CRST_THETA3_B);
	//	number f4t4 = _f4(t4, rnamodel.RNA_CRST_THETA4_T0, rnamodel.RNA_CRST_THETA4_TS, rnamodel.RNA_CRST_THETA4_TC, rnamodel.RNA_CRST_THETA4_A, rnamodel.RNA_CRST_THETA4_B) +
		//		_f4(PI - t4, rnamodel.RNA_CRST_THETA4_T0, rnamodel.RNA_CRST_THETA4_TS, rnamodel.RNA_CRST_THETA4_TC, rnamodel.RNA_CRST_THETA4_A, rnamodel.RNA_CRST_THETA4_B);
		number f4t7 = _f4(t7, rnamodel.RNA_CRST_THETA7_T0, rnamodel.RNA_CRST_THETA7_TS, rnamodel.RNA_CRST_THETA7_TC, rnamodel.RNA_CRST_THETA7_A, rnamodel.RNA_CRST_THETA7_B) +
				_f4(PI - t7, rnamodel.RNA_CRST_THETA7_T0, rnamodel.RNA_CRST_THETA7_TS, rnamodel.RNA_CRST_THETA7_TC, rnamodel.RNA_CRST_THETA7_A, rnamodel.RNA_CRST_THETA7_B);
		number f4t8 = _f4(t8, rnamodel.RNA_CRST_THETA8_T0, rnamodel.RNA_CRST_THETA8_TS, rnamodel.RNA_CRST_THETA8_TC, rnamodel.RNA_CRST_THETA8_A, rnamodel.RNA_CRST_THETA8_B) +
				_f4(PI - t8, rnamodel.RNA_CRST_THETA8_T0, rnamodel.RNA_CRST_THETA8_TS, rnamodel.RNA_CRST_THETA8_TC, rnamodel.RNA_CRST_THETA8_A, rnamodel.RNA_CRST_THETA8_B);

		number cstk_energy = f2 * f4t1 * f4t2 * f4t3 *  f4t7 * f4t8;

		if(cstk_energy < (number) 0) {
			// derivatives called at the relevant arguments
			number f2D = _f2D(rcstackmod, CRST_F2);
			number f4t1Dsin = -_f4Dsin(t1, rnamodel.RNA_CRST_THETA1_T0, rnamodel.RNA_CRST_THETA1_TS, rnamodel.RNA_CRST_THETA1_TC, rnamodel.RNA_CRST_THETA1_A, rnamodel.RNA_CRST_THETA1_B);
			number f4t2Dsin = -_f4Dsin(t2, rnamodel.RNA_CRST_THETA2_T0, rnamodel.RNA_CRST_THETA2_TS, rnamodel.RNA_CRST_THETA2_TC, rnamodel.RNA_CRST_THETA2_A, rnamodel.RNA_CRST_THETA2_B);
			number f4t3Dsin = _f4Dsin(t3, rnamodel.RNA_CRST_THETA3_T0, rnamodel.RNA_CRST_THETA3_TS, rnamodel.RNA_CRST_THETA3_TC, rnamodel.RNA_CRST_THETA3_A, rnamodel.RNA_CRST_THETA3_B);
			//number f4t4Dsin = _f4Dsin(t4, rnamodel.RNA_CRST_THETA4_T0, rnamodel.RNA_CRST_THETA4_TS, rnamodel.RNA_CRST_THETA4_TC, rnamodel.RNA_CRST_THETA4_A, rnamodel.RNA_CRST_THETA4_B) -
			//		_f4Dsin(PI - t4, rnamodel.RNA_CRST_THETA4_T0, rnamodel.RNA_CRST_THETA4_TS, rnamodel.RNA_CRST_THETA4_TC, rnamodel.RNA_CRST_THETA4_A, rnamodel.RNA_CRST_THETA4_B);
			number f4t7Dsin = -_f4Dsin(t7, rnamodel.RNA_CRST_THETA7_T0, rnamodel.RNA_CRST_THETA7_TS, rnamodel.RNA_CRST_THETA7_TC, rnamodel.RNA_CRST_THETA7_A, rnamodel.RNA_CRST_THETA7_B) +
					_f4Dsin(PI - t7, rnamodel.RNA_CRST_THETA7_T0, rnamodel.RNA_CRST_THETA7_TS, rnamodel.RNA_CRST_THETA7_TC, rnamodel.RNA_CRST_THETA7_A, rnamodel.RNA_CRST_THETA7_B);
			number f4t8Dsin = _f4Dsin(t8, rnamodel.RNA_CRST_THETA8_T0, rnamodel.RNA_CRST_THETA8_TS, rnamodel.RNA_CRST_THETA8_TC, rnamodel.RNA_CRST_THETA8_A, rnamodel.RNA_CRST_THETA8_B) -
					_f4Dsin(PI - t8, rnamodel.RNA_CRST_THETA8_T0, rnamodel.RNA_CRST_THETA8_TS, rnamodel.RNA_CRST_THETA8_TC, rnamodel.RNA_CRST_THETA8_A, rnamodel.RNA_CRST_THETA8_B);

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
			//Ttmp -= _cross<number, number4>(a3, b3) * (-cstk_energy * f4t4Dsin / f4t4);

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
	
	number4 rstack = r + qpos_stack - ppos_stack;
	number rstackmodsqr = CUDA_DOT(rstack, rstack);
	if(SQR(rnamodel.RNA_CXST_RCLOW) < rstackmodsqr && rstackmodsqr < SQR(rnamodel.RNA_CXST_RCHIGH)) {
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
		//number4 rbackboneref = r + rnamodel.RNA_POS_BACK * b1 - rnamodel.RNA_POS_BACK * a1;
		number rbackmod = _module<number, number4>(rbackbone);
		number4 rbackbonedir = rbackbone / rbackmod;

		number cosphi3 = CUDA_DOT(rstackdir, (_cross<number, number4>(rbackbonedir, a1)));
		number cosphi4 = CUDA_DOT(rstackdir, (_cross<number, number4>(rbackbonedir, b1)));

	 	// functions called at their relevant arguments
		number f2 = _f2(rstackmod, CXST_F2);
		number f4t1 = _f4(t1, rnamodel.RNA_CXST_THETA1_T0, rnamodel.RNA_CXST_THETA1_TS, rnamodel.RNA_CXST_THETA1_TC, rnamodel.RNA_CXST_THETA1_A, rnamodel.RNA_CXST_THETA1_B) +
				_f4(2 * PI - t1, rnamodel.RNA_CXST_THETA1_T0, rnamodel.RNA_CXST_THETA1_TS, rnamodel.RNA_CXST_THETA1_TC, rnamodel.RNA_CXST_THETA1_A, rnamodel.RNA_CXST_THETA1_B);
		number f4t4 = _f4(t4, rnamodel.RNA_CXST_THETA4_T0, rnamodel.RNA_CXST_THETA4_TS, rnamodel.RNA_CXST_THETA4_TC, rnamodel.RNA_CXST_THETA4_A, rnamodel.RNA_CXST_THETA4_B);
		number f4t5 = _f4(t5, rnamodel.RNA_CXST_THETA5_T0, rnamodel.RNA_CXST_THETA5_TS, rnamodel.RNA_CXST_THETA5_TC, rnamodel.RNA_CXST_THETA5_A, rnamodel.RNA_CXST_THETA5_B) +
				_f4(PI - t5, rnamodel.RNA_CXST_THETA5_T0, rnamodel.RNA_CXST_THETA5_TS, rnamodel.RNA_CXST_THETA5_TC, rnamodel.RNA_CXST_THETA5_A, rnamodel.RNA_CXST_THETA5_B);
		number f4t6 = _f4(t6, rnamodel.RNA_CXST_THETA6_T0, rnamodel.RNA_CXST_THETA6_TS, rnamodel.RNA_CXST_THETA6_TC, rnamodel.RNA_CXST_THETA6_A, rnamodel.RNA_CXST_THETA6_B) +
				_f4(PI - t6, rnamodel.RNA_CXST_THETA6_T0, rnamodel.RNA_CXST_THETA6_TS, rnamodel.RNA_CXST_THETA6_TC, rnamodel.RNA_CXST_THETA6_A, rnamodel.RNA_CXST_THETA6_B);
		number f5cosphi3 = _f5(cosphi3, RNA_CXST_F5_PHI3);

		number f5cosphi4 = _f5(cosphi4, RNA_CXST_F5_PHI4);


		number cxst_energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5cosphi4;

		//if(grooving)  printf("AADEVICE: Energy =  %f,  f2 = %f, f4t1=%f, f4t4=%f, f4t5=%f, f4t6=%f, f5cosphi3=%f, f5cosphi4=%f \n",cxst_energy,f2 , f4t1 , f4t4 ,f4t5 , f4t6 ,f5cosphi3 , f5cosphi4);

		//
		if(cxst_energy < (number) 0) {
			// derivatives called at the relevant arguments
			number f2D = _f2D(rstackmod, RNA_CXST_F2);
			number f4t1Dsin = -_f4Dsin(t1, rnamodel.RNA_CXST_THETA1_T0, rnamodel.RNA_CXST_THETA1_TS, rnamodel.RNA_CXST_THETA1_TC, rnamodel.RNA_CXST_THETA1_A, rnamodel.RNA_CXST_THETA1_B) +
					_f4Dsin(2 * PI - t1, rnamodel.RNA_CXST_THETA1_T0, rnamodel.RNA_CXST_THETA1_TS, rnamodel.RNA_CXST_THETA1_TC, rnamodel.RNA_CXST_THETA1_A, rnamodel.RNA_CXST_THETA1_B);
			number f4t4Dsin =  _f4Dsin(t4, rnamodel.RNA_CXST_THETA4_T0, rnamodel.RNA_CXST_THETA4_TS, rnamodel.RNA_CXST_THETA4_TC, rnamodel.RNA_CXST_THETA4_A, rnamodel.RNA_CXST_THETA4_B);
			number f4t5Dsin =  _f4Dsin(t5, rnamodel.RNA_CXST_THETA5_T0, rnamodel.RNA_CXST_THETA5_TS, rnamodel.RNA_CXST_THETA5_TC, rnamodel.RNA_CXST_THETA5_A, rnamodel.RNA_CXST_THETA5_B) -
					_f4Dsin(PI - t5, rnamodel.RNA_CXST_THETA5_T0, rnamodel.RNA_CXST_THETA5_TS, rnamodel.RNA_CXST_THETA5_TC, rnamodel.RNA_CXST_THETA5_A, rnamodel.RNA_CXST_THETA5_B);
			number f4t6Dsin = -_f4Dsin(t6, rnamodel.RNA_CXST_THETA6_T0, rnamodel.RNA_CXST_THETA6_TS, rnamodel.RNA_CXST_THETA6_TC, rnamodel.RNA_CXST_THETA6_A, rnamodel.RNA_CXST_THETA6_B) +
					_f4Dsin(PI - t6, rnamodel.RNA_CXST_THETA6_T0, rnamodel.RNA_CXST_THETA6_TS, rnamodel.RNA_CXST_THETA6_TC, rnamodel.RNA_CXST_THETA6_A, rnamodel.RNA_CXST_THETA6_B);
			number f5Dcosphi3 = _f5D(cosphi3, RNA_CXST_F5_PHI3);
			number f5Dcosphi4 = _f5D(cosphi4, RNA_CXST_F5_PHI4);

			// RADIAL PART
			Ftmp = -rstackdir * (cxst_energy * f2D / f2);

			//if(grooving)  printf("DEVICE radial force: Ftmp = (%f,%f,%f) with stackdir = (%f,%f,%f), a3= (%f,%f,%f) \n",Ftmp.x,Ftmp.y,Ftmp.z,rstackdir.x,rstackdir.y,rstackdir.z,a3.x,a3.y,a3.z);
			//if(grooving)  printf("AADEVICE:, with rstackmod=%g, Force terms f2 = %g, f4t1=%g, f4t4=%g, f4t5=%g, f4t6=%g, f5cosphi3=%g, f5cosphi4=%g \n",rstackmod,f2D , f4t1Dsin , f4t4Dsin ,f4t5Dsin , f4t6Dsin ,f5Dcosphi3 , f5Dcosphi4);

			// THETA1; t1 = LRACOS (-a1 * b1);
			Ttmp -= _cross<number, number4>(a1, b1) * (-cxst_energy * f4t1Dsin / f4t1);

			// TETA4; t4 = LRACOS (a3 * b3);
			Ttmp -= _cross<number, number4>(a3, b3) * (-cxst_energy * f4t4Dsin / f4t4);

			// THETA5; t5 = LRACOS ( a3 * rstackdir);
			number part = cxst_energy * f4t5Dsin / f4t5;
			Ftmp += (a3 - rstackdir * cosf(t5)) / rstackmod * part;
			Ttmp -= _cross<number, number4>(rstackdir, a3) * part;
			//if(grooving)  printf("DEVICE theta5 force: Ftmp = (%f,%f,%f) with stackdir = (%f,%f,%f), a3= (%f,%f,%f) \n",Ftmp.x,Ftmp.y,Ftmp.z,rstackdir.x,rstackdir.y,rstackdir.z,a3.x,a3.y,a3.z);


			// THETA6; t6 = LRACOS (-b3 * rstackdir);
			Ftmp += (b3 + rstackdir * cosf(t6)) * (cxst_energy * f4t6Dsin / (f4t6 * rstackmod));


			//if(grooving)  printf("Before: Ttmp = (%f,%f,%f)\n",Ttmp.x,Ttmp.y,Ttmp.z);

			Ttmp -= _cross<number, number4>(ppos_stack, Ftmp);
			//if(grooving)  printf("AADEVICE: Ttmp = (%f,%f,%f)\n",Ttmp.x,Ttmp.y,Ttmp.z);


			number4 rba1 =  _cross<number,number4>(rbackbonedir,a1);
			number4 a1rs =  _cross<number,number4>(a1,rstackdir);
			number rb_dot_a1rs = CUDA_DOT(rbackbonedir , a1rs);
			number rs_dot_rba1 = CUDA_DOT(rstackdir , rba1);

			number force_c = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi4 * f5Dcosphi3;

			number4 forcestack = -(rba1  - rstackdir * rs_dot_rba1) *  (force_c / rstackmod);
			number4 forceback = -(a1rs - rbackbonedir * rb_dot_a1rs) * (force_c   / rbackmod);

			number4 myforce =  forcestack;
			myforce +=  forceback;


			// for p
			number4 mytorque1  =-  _cross<number,number4>(ppos_stack, forcestack ); //a1 * (p->int_centers[RNANucleotide<number>::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide<number>::STACK] * a1 )
			mytorque1  -=  _cross<number,number4>(ppos_back, forceback) ;  // p->int_centers[RNANucleotide<number>::BACK].cross(forceback);

			// for q
			number4 mytorque2  = _cross<number,number4>(qpos_stack,forcestack); //  q->int_centers[RNANucleotide<number>::STACK].cross( forcestack ); //a1 * (p->int_centers[RNANucleotide<number>::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide<number>::STACK] * a1 )
			mytorque2 +=   _cross<number,number4>(qpos_back,forceback); // q->int_centers[RNANucleotide<number>::BACK].cross(forceback);

			mytorque1 -=  force_c  *   _cross<number,number4>(a1,  _cross<number,number4>(rstackdir,rbackbonedir)); // a1.cross(rstackdir.cross(rbackbonedir));  // rstackdir * (rbackbonedir * a1) - rbackbonedir * (rstackdir * a1);


			force_c = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5Dcosphi4;
		    rba1 =     _cross<number,number4>(rbackbonedir,b1); // rbackbonedir.cross(b1);
			a1rs =     _cross<number,number4>(b1,rstackdir);   //    b1.cross(rstackdir);
			rb_dot_a1rs = CUDA_DOT(rbackbonedir , a1rs);
			rs_dot_rba1 = CUDA_DOT(rstackdir , rba1);

			forcestack = -(rba1  - rstackdir * rs_dot_rba1) *  (force_c / rstackmod);
		    forceback = -(a1rs - rbackbonedir * rb_dot_a1rs) * (force_c   / rbackmod);

			myforce +=  forcestack;
			myforce +=  forceback;

									// for p
			mytorque2  +=   _cross<number,number4>(qpos_stack,forcestack); // q->int_centers[RNANucleotide<number>::STACK].cross( forcestack ); //a1 * (p->int_centers[RNANucleotide<number>::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide<number>::STACK] * a1 )
			mytorque2  +=  _cross<number,number4>(qpos_back,forceback); //  q->int_centers[RNANucleotide<number>::BACK].cross(forceback);

									// for q
			mytorque1  -=  _cross<number,number4>(ppos_stack,forcestack); // p->int_centers[RNANucleotide<number>::STACK].cross( forcestack ); //a1 * (p->int_centers[RNANucleotide<number>::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide<number>::STACK] * a1 )
			mytorque1 -=  _cross<number,number4>(ppos_back,forceback); //p->int_centers[RNANucleotide<number>::BACK].cross(forceback);

			mytorque2 -=  force_c  *  _cross<number,number4>(b1,rbackbonedir); // b1.cross(rstackdir.cross(rbackbonedir));  // rstackdir * (rbackbonedir * a1) - rbackbonedir * (rstackdir * a1);



			Ftmp -= myforce;
			Ttmp += mytorque1;



			//if(grooving)  printf("XXX DEVICE: Ftmp=( %f %f %f ), Ttmp=( %f %f %f ), Energy =  %f,  f2 = %f, f4t1=%f, f4t4=%f, f4t5=%f, f4t6=%f, f5cosphi3=%f, f5cosphi4=%f \n",Ftmp.x,Ftmp.y,Ftmp.z, Ttmp.x, Ttmp.y, Ttmp.z,cxst_energy,f2 , f4t1 , f4t4 ,f4t5 , f4t6 ,f5cosphi3 , f5cosphi4);

			//if(grooving)  printf("XXX DEVICE: Ftmp=( %f %f %f ), Ttmp=( %f %f %f )\n",Ftmp.x,Ftmp.y,Ftmp.z, Ttmp.x, Ttmp.y, Ttmp.z);

			//if(grooving && (f5Dcosphi4 != 0 || f5Dcosphi3 != 0))  printf("DEVICE: f5Dcosphi4: %f f5Dcosphi3: %f  \n",f5Dcosphi4,f5Dcosphi3);

			//if(grooving)  printf("Energy =  %f,  f2 = %f, f4t1=%f, f4t4=%f, f4t5=%f, f4t6=%f, f5cosphi3=%f, f5cosphi4=%f \n",Ftmp.x,Ftmp.y,Ftmp.z, Ttmp.x, Ttmp.y, Ttmp.z,cxst_energy,f2 , f4t1 , f4t4 ,f4t5 , f4t6 ,f5cosphi3 , f5cosphi4);



//
//
//			// COSPHI3
//			number rbackrefmodcub = rbackrefmod * rbackrefmod * rbackrefmod;
//
//			//number a1b1 = a1 * b1;
//			number a2b1 = CUDA_DOT(a2, b1);
//			number a3b1 = CUDA_DOT(a3, b1);
//			number ra1 = CUDA_DOT(rstackdir, a1);
//			number ra2 = CUDA_DOT(rstackdir, a2);
//			number ra3 = CUDA_DOT(rstackdir, a3);
//			number rb1 = CUDA_DOT(rstackdir, b1);
//
//			number parentesi = (ra3 * a2b1 - ra2 * a3b1);
//			number dcdr    = -rnamodel.RNA_GAMMA * parentesi * (rnamodel.RNA_GAMMA * (ra1 - rb1) + rstackmod) / rbackrefmodcub;
//			number dcda1b1 =  rnamodel.RNA_GAMMA * SQR(rnamodel.RNA_GAMMA) * parentesi / rbackrefmodcub;
//			number dcda2b1 =  rnamodel.RNA_GAMMA * ra3 / rbackrefmod;
//			number dcda3b1 = -rnamodel.RNA_GAMMA * ra2 / rbackrefmod;
//			number dcdra1  = -SQR(rnamodel.RNA_GAMMA) * parentesi * rstackmod / rbackrefmodcub;
//			number dcdra2  = -rnamodel.RNA_GAMMA * a3b1 / rbackrefmod;
//			number dcdra3  =  rnamodel.RNA_GAMMA * a2b1 / rbackrefmod;
//			number dcdrb1  = -dcdra1;
//
//			part = cxst_energy * 2 * f5cosphi3D / f5cosphi3;
//
//			Ftmp -= part * (rstackdir * dcdr +
//						    ((a1 - rstackdir * ra1) * dcdra1 +
//							(a2 - rstackdir * ra2) * dcdra2 +
//							(a3 - rstackdir * ra3) * dcdra3 +
//							(b1 - rstackdir * rb1) * dcdrb1) / rstackmod);
//
//			Ttmp += part * (_cross<number, number4>(rstackdir, a1) * dcdra1 +
//						    _cross<number, number4>(rstackdir, a2) * dcdra2 +
//						    _cross<number, number4>(rstackdir ,a3) * dcdra3);
//
//			Ttmp -= part * (_cross<number, number4>(a1, b1) * dcda1b1 +
//						    _cross<number, number4>(a2, b1) * dcda2b1 +
//						    _cross<number, number4>(a3, b1) * dcda3b1);
//



			Ftmp.w = cxst_energy;
			F -= Ftmp;



		}
	}

	//DEBYE-HUCKEL
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

	//if(grooving)  printf("XXX DEVICE all pairs: Ftmp=( %f %f %f ), Ttmp=( %f %f %f ), Energy =  %f \n",Ftmp.x,Ftmp.y,Ftmp.z, Ttmp.x, Ttmp.y, Ttmp.z,Ftmp.w);


	// this component stores the energy due to hydrogen bonding
	T.w = old_Tw + hb_energy;
}

// forces + second step without lists
template <typename number, typename number4>
__global__ void rna_forces(number4 *poss, GPU_quat<number> *orientations, number4 *forces, number4 *torques, LR_bonds *bonds, bool average, bool use_debye_huckel, bool mismatch_repulsion) {
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
						    qpos, b1, b2, b3, F, T, average);


	}

	if(bs.n5 != P_INVALID) {
		number4 qpos = poss[bs.n5];

		number4 b1, b2, b3;
		get_vectors_from_quat<number,number4>(orientations[bs.n5], b1, b2, b3);
		_bonded_part<number, number4, false>(qpos, b1, b2, b3,
						     ppos, a1, a2, a3, F, T, average);

	}

	const int type = get_particle_type<number, number4>(ppos);
	T.w = (number) 0;
	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND && bs.n3 != j && bs.n5 != j) {
			const number4 qpos = poss[j];
			number4 b1, b2, b3;
			get_vectors_from_quat<number,number4>(orientations[j], b1, b2, b3);
			LR_bonds qbonds = bonds[j];

			_particle_particle_interaction<number, number4>(ppos, a1, a2, a3, qpos, b1, b2, b3, F, T, average,use_debye_huckel, mismatch_repulsion, bs, qbonds);


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
__global__ void rna_forces_edge_nonbonded(number4 *poss, GPU_quat<number> *orientations, number4 *forces, number4 *torques, edge_bond *edge_list, int n_edges,LR_bonds *bonds, bool average,bool use_debye_huckel, bool mismatch_repulsion) {
	if(IND >= n_edges) return;

	number4 dF = make_number4<number, number4>(0, 0, 0, 0);
	number4 dT = make_number4<number, number4>(0, 0, 0, 0);

	edge_bond b = edge_list[IND];

	// get info for particle 1
	number4 ppos = poss[b.from];
	number4 a1, a2, a3;
	get_vectors_from_quat<number,number4>(orientations[b.from], a1, a2, a3);

	// get info for particle 2
	number4 qpos = poss[b.to];
	number4 b1, b2, b3;
	get_vectors_from_quat<number,number4>(orientations[b.to], b1, b2, b3);

	LR_bonds pbonds = bonds[b.from];
	LR_bonds qbonds = bonds[b.to];
	_particle_particle_interaction<number, number4>(ppos, a1, a2, a3, qpos, b1, b2, b3, dF, dT, average, use_debye_huckel, mismatch_repulsion , pbonds, qbonds);

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
__global__ void rna_forces_edge_bonded(number4 *poss, GPU_quat<number> *orientations,  number4 *forces, number4 *torques, LR_bonds *bonds, bool average) {
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
	number4 a1, a2, a3;
	get_vectors_from_quat<number,number4>(orientations[IND], a1, a2, a3);

	if(bs.n3 != P_INVALID) {
		number4 qpos = poss[bs.n3];
		number4 b1,b2,b3;
		get_vectors_from_quat<number,number4>(orientations[bs.n3], b1, b2, b3);


		//grooving =false;
		_bonded_part<number, number4, true>(ppos, a1, a2, a3, qpos, b1, b2, b3, dF, dT, average);
	}
	if(bs.n5 != P_INVALID) {
		number4 qpos = poss[bs.n5];
		number4 b1,b2,b3;
		get_vectors_from_quat<number,number4>(orientations[bs.n5], b1, b2, b3);



		_bonded_part<number, number4, false>(qpos, b1, b2, b3, ppos, a1, a2, a3, dF, dT, average);


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
__global__ void rna_forces(number4 *poss, GPU_quat<number> *orientations,  number4 *forces, number4 *torques, int *matrix_neighs, int *number_neighs, LR_bonds *bonds, bool average, bool use_debye_huckel, bool mismatch_repulsion) {
	if(IND >= MD_N[0]) return;

	//number4 F = make_number4<number, number4>(0, 0, 0, 0);
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
						    qpos, b1, b2, b3, F, T, average);


	}
	if(bs.n5 != P_INVALID) {
		number4 qpos = poss[bs.n5];
		number4 b1,b2,b3;
		get_vectors_from_quat<number,number4>(orientations[bs.n5], b1, b2, b3);


		_bonded_part<number, number4, false>(qpos, b1, b2, b3,
						     ppos, a1, a2, a3, F, T, average);
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
		_particle_particle_interaction<number, number4>(ppos, a1, a2, a3, qpos, b1, b2, b3, F, T, average,use_debye_huckel, mismatch_repulsion, pbonds, qbonds);
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
__global__ void rna_hb_op_precalc(number4 *poss, GPU_quat<number> *orientations, int *op_pairs1, int *op_pairs2, float *hb_energies, int n_threads, bool *region_is_nearhb)
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

	number4 ppos_base = rnamodel.RNA_POS_BASE * a1;
	number4 qpos_base = rnamodel.RNA_POS_BASE * b1;

	// HYDROGEN BONDING

	/*
	if(!average) //allow for wobble bp
	{
			if( int_type == 4 && ( (qtype == N_T && ptype == N_G ) || (qtype == N_G && ptype == N_T)  ))
				is_pair = 1;
	}*/

	int is_pair = int(int_type == 3);
	number hb_energy = (number) 0;
	number4 rhydro = r + qpos_base - ppos_base;
	number rhydromodsqr = CUDA_DOT(rhydro, rhydro);
	if(is_pair && SQR(rnamodel.RNA_HYDR_RCLOW) < rhydromodsqr && rhydromodsqr < SQR(rnamodel.RNA_HYDR_RCHIGH)) {
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
		number f1 = hb_multi * _f1(rhydromod, RNA_HYDR_F1, ptype, qtype);
		number f4t1 = _f4(t1, rnamodel.RNA_HYDR_THETA1_T0, rnamodel.RNA_HYDR_THETA1_TS, rnamodel.RNA_HYDR_THETA1_TC, rnamodel.RNA_HYDR_THETA1_A, rnamodel.RNA_HYDR_THETA1_B);
		number f4t2 = _f4(t2, rnamodel.RNA_HYDR_THETA2_T0, rnamodel.RNA_HYDR_THETA2_TS, rnamodel.RNA_HYDR_THETA2_TC, rnamodel.RNA_HYDR_THETA2_A, rnamodel.RNA_HYDR_THETA2_B);
		number f4t3 = _f4(t3, rnamodel.RNA_HYDR_THETA3_T0, rnamodel.RNA_HYDR_THETA3_TS, rnamodel.RNA_HYDR_THETA3_TC, rnamodel.RNA_HYDR_THETA3_A, rnamodel.RNA_HYDR_THETA3_B);
		number f4t4 = _f4(t4, rnamodel.RNA_HYDR_THETA4_T0, rnamodel.RNA_HYDR_THETA4_TS, rnamodel.RNA_HYDR_THETA4_TC, rnamodel.RNA_HYDR_THETA4_A, rnamodel.RNA_HYDR_THETA4_B);
		number f4t7 = _f4(t7, rnamodel.RNA_HYDR_THETA7_T0, rnamodel.RNA_HYDR_THETA7_TS, rnamodel.RNA_HYDR_THETA7_TC, rnamodel.RNA_HYDR_THETA7_A, rnamodel.RNA_HYDR_THETA7_B);
		number f4t8 = _f4(t8, rnamodel.RNA_HYDR_THETA8_T0, rnamodel.RNA_HYDR_THETA8_TS, rnamodel.RNA_HYDR_THETA8_TC, rnamodel.RNA_HYDR_THETA8_A, rnamodel.RNA_HYDR_THETA8_B);

		hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
	}
	//END HYDROGEN BONDING

	hb_energies[IND] = hb_energy;
}

template <typename number, typename number4>
__global__ void rna_near_hb_op_precalc(number4 *poss, GPU_quat<number> *orientations, int *op_pairs1, int *op_pairs2, bool *nearly_bonded_array, int n_threads, bool *region_is_nearhb)
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

	number4 ppos_base = rnamodel.RNA_POS_BASE * a1;
	number4 qpos_base = rnamodel.RNA_POS_BASE * b1;

	// HYDROGEN BONDING

	int is_pair = int(int_type == 3);
	/*if(!average) //allow for wobble bp
	{
			if( int_type == 4 && ( (qtype == N_T && ptype == N_G ) || (qtype == N_G && ptype == N_T)  ))
				is_pair = 1;
	}*/

	number4 rhydro = r + qpos_base - ppos_base;
	number rhydromodsqr = CUDA_DOT(rhydro, rhydro);

	int total_nonzero;
	bool nearly_bonded = false;

	if(is_pair && SQR(rnamodel.RNA_HYDR_RCLOW) < rhydromodsqr && rhydromodsqr < SQR(rnamodel.RNA_HYDR_RCHIGH)) {
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
		number f1 = hb_multi * _f1(rhydromod, RNA_HYDR_F1, ptype, qtype);
		number f4t1 = _f4(t1, rnamodel.RNA_HYDR_THETA1_T0, rnamodel.RNA_HYDR_THETA1_TS, rnamodel.RNA_HYDR_THETA1_TC, rnamodel.RNA_HYDR_THETA1_A, rnamodel.RNA_HYDR_THETA1_B);
		number f4t2 = _f4(t2, rnamodel.RNA_HYDR_THETA2_T0, rnamodel.RNA_HYDR_THETA2_TS, rnamodel.RNA_HYDR_THETA2_TC, rnamodel.RNA_HYDR_THETA2_A, rnamodel.RNA_HYDR_THETA2_B);
		number f4t3 = _f4(t3, rnamodel.RNA_HYDR_THETA3_T0, rnamodel.RNA_HYDR_THETA3_TS, rnamodel.RNA_HYDR_THETA3_TC, rnamodel.RNA_HYDR_THETA3_A, rnamodel.RNA_HYDR_THETA3_B);
		number f4t4 = _f4(t4, rnamodel.RNA_HYDR_THETA4_T0, rnamodel.RNA_HYDR_THETA4_TS, rnamodel.RNA_HYDR_THETA4_TC, rnamodel.RNA_HYDR_THETA4_A, rnamodel.RNA_HYDR_THETA4_B);
		number f4t7 = _f4(t7, rnamodel.RNA_HYDR_THETA7_T0, rnamodel.RNA_HYDR_THETA7_TS, rnamodel.RNA_HYDR_THETA7_TC, rnamodel.RNA_HYDR_THETA7_A, rnamodel.RNA_HYDR_THETA7_B);
		number f4t8 = _f4(t8, rnamodel.RNA_HYDR_THETA8_T0, rnamodel.RNA_HYDR_THETA8_TS, rnamodel.RNA_HYDR_THETA8_TC, rnamodel.RNA_HYDR_THETA8_A, rnamodel.RNA_HYDR_THETA8_B);

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
__global__ void rna_dist_op_precalc(number4 *poss, GPU_quat<number> *orientations, int *op_pairs1, int *op_pairs2, number *op_dists, int n_threads)
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

	number4 ppos_base = rnamodel.RNA_POS_BASE * a1;
	number4 qpos_base = rnamodel.RNA_POS_BASE * b1;
	
	number4 rbase = r + qpos_base - ppos_base;
	op_dists[IND] = _module<number, number4>(rbase);
}
