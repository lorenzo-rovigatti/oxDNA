__constant__ float MD_sqr_verlet_skin[1];
/* System constants */
__constant__ int MD_N[1];
__constant__ float MD_dt[1];

#include "../cuda_utils/CUDA_lr_common.cuh"

// this function modifies L
__device__ LR_GPU_matrix<double> _get_updated_orientation(LR_double4 &L, const LR_GPU_matrix<double> &old_o) {
	double norm = sqrt(CUDA_DOT(L, L));
	L.x /= norm;
	L.y /= norm;
	L.z /= norm;

	double sintheta, costheta;
	sincos(MD_dt[0] * norm, &sintheta, &costheta);

	double olcos = 1. - costheta;
	double xyo = L.x * L.y * olcos;
	double xzo = L.x * L.z * olcos;
	double yzo = L.y * L.z * olcos;
	double xsin = L.x * sintheta;
	double ysin = L.y * sintheta;
	double zsin = L.z * sintheta;

	LR_GPU_matrix<double> R = {SQR(L.x) * olcos + costheta, xyo - zsin, xzo + ysin,
						  xyo + zsin, SQR(L.y) * olcos + costheta, yzo - xsin,
						  xzo - ysin, yzo + xsin, SQR(L.z) * olcos + costheta};

	return _matrix_matrix_product<double>(old_o, R);
}

__global__ void first_step_mixed(float4 *poss, LR_GPU_matrix<float> *orientations, LR_double4 *possd, LR_GPU_matrix<double> *orientationsd, float4 *list_poss, LR_double4 *velsd, LR_double4 *Lsd, float4 *forces, float4 *torques, bool *are_lists_old, bool any_rigid_body) {
	if(IND >= MD_N[0]) return;

	float4 F = forces[IND];

	LR_double4 v = velsd[IND];

	v.x += F.x * MD_dt[0] * 0.5f;
	v.y += F.y * MD_dt[0] * 0.5f;
	v.z += F.z * MD_dt[0] * 0.5f;

	velsd[IND] = v;

	LR_double4 r = possd[IND];

	r.x += v.x * MD_dt[0];
	r.y += v.y * MD_dt[0];
	r.z += v.z * MD_dt[0];

	possd[IND] = r;

	float4 rf = make_float4(r.x, r.y, r.z, r.w);
	poss[IND] = rf;

	if(any_rigid_body) {
		float4 T = torques[IND];
		LR_double4 L = Lsd[IND];
		
		L.x += T.x * MD_dt[0] * 0.5;
		L.y += T.y * MD_dt[0] * 0.5;
		L.z += T.z * MD_dt[0] * 0.5;
		
		Lsd[IND] = L;
		
		const LR_GPU_matrix<double> new_o = _get_updated_orientation(L, orientationsd[IND]);
		orientationsd[IND] = new_o;
		
		const LR_GPU_matrix<float> new_of = {(float)new_o.e[0], (float)new_o.e[1], (float)new_o.e[2],
									  (float)new_o.e[3], (float)new_o.e[4], (float)new_o.e[5],
									  (float)new_o.e[6], (float)new_o.e[7], (float)new_o.e[8]};
		orientations[IND] = new_of;
	}

	// do verlet lists need to be updated?
	if(quad_distance<float, float4>(rf, list_poss[IND]) > MD_sqr_verlet_skin[0]) are_lists_old[0] = true;
}

__global__ void second_step_mixed(LR_double4 *velsd, LR_double4 *Lsd, float4 *forces, float4 *torques, bool any_rigid_body) {
	if(IND >= MD_N[0]) return;

	float4 F = forces[IND];
	LR_double4 v = velsd[IND];

	v.x += (F.x * MD_dt[0] * 0.5f);
	v.y += (F.y * MD_dt[0] * 0.5f);
	v.z += (F.z * MD_dt[0] * 0.5f);

	velsd[IND] = v;

	if(any_rigid_body) {
		float4 T = torques[IND];
		LR_double4 L = Lsd[IND];
		
		L.x += (T.x * MD_dt[0] * 0.5f);
		L.y += (T.y * MD_dt[0] * 0.5f);
		L.z += (T.z * MD_dt[0] * 0.5f);
		
		Lsd[IND] = L;
	}
}

__global__ void float4_to_LR_double4(float4 *src, LR_double4 *dest) {
	if(IND >= MD_N[0]) return;

	float4 tmp = src[IND];
	dest[IND] = make_LR_double4(tmp);
}

__global__ void LR_double4_to_float4(LR_double4 *src, float4 *dest) {
	if(IND >= MD_N[0]) return;

	LR_double4 tmp = src[IND];
	dest[IND] = make_float4((float)tmp.x, (float)tmp.y, (float)tmp.z, tmp.w);
}

__global__ void float4_to_LR_double4(LR_GPU_matrix<float> *src, LR_GPU_matrix<double> *dest) {
	if(IND >= MD_N[0]) return;

	LR_GPU_matrix<float> tmp = src[IND];
	LR_GPU_matrix<double> res = {(double)tmp.e[0], (double)tmp.e[1], (double)tmp.e[2],
				(double)tmp.e[3], (double)tmp.e[4], (double)tmp.e[5],
				(double)tmp.e[6], (double)tmp.e[7], (double)tmp.e[8]};
	dest[IND] = res;
}

__global__ void LR_double4_to_float4(LR_GPU_matrix<double> *src, LR_GPU_matrix<float> *dest) {
	if(IND >= MD_N[0]) return;

	LR_GPU_matrix<double> tmp = src[IND];
	LR_GPU_matrix<float> res = {(float)tmp.e[0], (float)tmp.e[1], (float)tmp.e[2],
				(float)tmp.e[3], (float)tmp.e[4], (float)tmp.e[5],
				(float)tmp.e[6], (float)tmp.e[7], (float)tmp.e[8]};
	dest[IND] = res;
}
