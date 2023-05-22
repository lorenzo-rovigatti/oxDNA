__constant__ float MD_sqr_verlet_skin[1];
__constant__ float MD_dt[1];
__constant__ int MD_N[1];

#include "../cuda_utils/CUDA_lr_common.cuh"

// this function modifies L
__device__ GPU_quat_double _get_updated_orientation(LR_double4 &L, GPU_quat_double &old_o) {
	double norm = sqrt(CUDA_DOT(L, L));
	L.x /= norm;
	L.y /= norm;
	L.z /= norm;
	
	double sintheta, costheta;
	sincos(MD_dt[0] * norm, &sintheta, &costheta);
	double qw = 0.5*sqrt(max(0., 2. + 2.*costheta));
	double winv = (double)1.0 /qw;
	GPU_quat_double R = {0.5*L.x*sintheta*winv, 0.5*L.y*sintheta*winv, 0.5*L.z*sintheta*winv, qw};
	
	return quat_multiply(old_o, R);
}

__global__ void first_step_mixed(float4 __restrict__ *poss, GPU_quat __restrict__ *orientations, LR_double4 __restrict__ *possd,
		GPU_quat_double __restrict__ *orientationsd, const float4 __restrict__ *list_poss, LR_double4 __restrict__ *velsd,
		LR_double4 __restrict__ *Lsd, const float4 __restrict__ *forces, const float4 __restrict__ *torques,
		bool *are_lists_old, bool any_rigid_body) {
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
		
		L.x += T.x * MD_dt[0] * 0.5f;
		L.y += T.y * MD_dt[0] * 0.5f;
		L.z += T.z * MD_dt[0] * 0.5f;
		
		Lsd[IND] = L;
		
		GPU_quat_double new_o = _get_updated_orientation(L, orientationsd[IND]);
		orientationsd[IND] = new_o;

		GPU_quat new_of = {(float)new_o.x, (float)new_o.y, (float)new_o.z, (float)new_o.w};
		orientations[IND] = new_of; 
	}

	// do verlet lists need to be updated?
	if(quad_distance(rf, list_poss[IND]) > MD_sqr_verlet_skin[0]) are_lists_old[0] = true;
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

__global__ void float4_to_LR_double4(GPU_quat *src, GPU_quat_double *dest) {
	if(IND >= MD_N[0]) return;

	GPU_quat tmp = src[IND];
	GPU_quat_double res = {(double)tmp.x, (double)tmp.y, (double)tmp.z, (double)tmp.w};
	dest[IND] = res;
}

__global__ void LR_double4_to_float4(GPU_quat_double *src, GPU_quat *dest) {
	if(IND >= MD_N[0]) return;

	GPU_quat_double tmp = src[IND];
	GPU_quat res = {(float)tmp.x, (float)tmp.y, (float)tmp.z, (float)tmp.w};
	dest[IND] = res;
}
