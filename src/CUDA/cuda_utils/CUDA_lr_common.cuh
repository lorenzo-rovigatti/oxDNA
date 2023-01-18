/*
 * CUDA_lr_common.cuh
 *
 *  Created on: Jun 3, 2011
 *      Author: rovigatti
 *
 *      I put these functions in a different file to tidy things up
 */

#ifndef CUDA_LR_COMMON
#define CUDA_LR_COMMON

#include "../CUDAUtils.h"

#include <curand_kernel.h>

/// Used to avoid syntax highlighting errors when coding with Eclipse
#ifdef __CDT_PARSER__
#define __global__
#define __device__
#define __shared__
#define IND 0
#endif

// in the case of pair-wise forces, which are counted twice, one should set half=true
// while three-body forces are counted only once, so half=false should be used.
template<bool half>
__device__ void _update_stress_tensor(CUDAStressTensor &st, const c_number4 &r, const c_number4 &force) {
	c_number factor = (half) ? 0.5f : 1.0f;

	st.e[0] -= r.x * force.x * factor;
	st.e[1] -= r.y * force.y * factor;
	st.e[2] -= r.z * force.z * factor;
	st.e[3] -= r.x * force.y * factor;
	st.e[4] -= r.x * force.z * factor;
	st.e[5] -= r.y * force.z * factor;
}

//This is the most commonly called quaternion to matrix conversion. 
__forceinline__ __device__ void get_vectors_from_quat(const GPU_quat &q, c_number4 &a1, c_number4 &a2, c_number4 &a3) {
	c_number sqx = q.x * q.x;
	c_number sqy = q.y * q.y;
	c_number sqz = q.z * q.z;
	c_number sqw = q.w * q.w;
	c_number xy = q.x * q.y;
	c_number xz = q.x * q.z;
	c_number xw = q.x * q.w;
	c_number yz = q.y * q.z;
	c_number yw = q.y * q.w;
	c_number zw = q.z * q.w;

	a1.x = (sqx - sqy - sqz + sqw);
	a2.x = 2 * (xy - zw);
	a3.x = 2 * (xz + yw);
	a1.y = 2 * (xy + zw);
	a2.y = (-sqx + sqy - sqz + sqw);
	a3.y = 2 * (yz - xw);
	a1.z = 2 * (xz - yw);
	a2.z = 2 * (yz + xw);
	a3.z = (-sqx - sqy + sqz + sqw);
}

template<typename t_quat>
__forceinline__ __device__ t_quat quat_multiply(const t_quat &a, const t_quat &b) {
	t_quat p;
	p.w = a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z;
	p.x = a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y;
	p.y = a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x;
	p.z = a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w;
	return p;
}

__forceinline__ __device__ void gaussian(curandState &state, float &outx, float &outy) {
	float r = sqrtf(-2.0f * logf(curand_uniform(&state)));
	float phi = 2.f * PI * curand_uniform(&state);

	outx = r * __cosf(phi);
	outy = r * __sinf(phi);
}

__forceinline__ __device__ void gaussian(curandState &state, double &outx, double &outy) {
	double r = sqrt(-2. * log(curand_uniform(&state)));
	double phi = 2 * M_PI * curand_uniform(&state);

	outx = r * cos(phi);
	outy = r * sin(phi);
}

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
/**
 * @brief found at { @link http://stackoverflow.com/questions/12626096/why-has-atomicadd-not-been-implemented-for-doubles }
 */
__forceinline__ __device__ double atomicAdd(double* address, double val) {
	unsigned long long int* address_as_ull = (unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;
	do {
		assumed = old;
		old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
	}
	while (assumed != old);

	return __longlong_as_double(old);
}
#endif

__forceinline__ __device__ void LR_atomicAdd(c_number4 *dst, const c_number4 &delta) {
	atomicAdd(&(dst->x), delta.x);
	atomicAdd(&(dst->y), delta.y);
	atomicAdd(&(dst->z), delta.z);
	atomicAdd(&(dst->w), delta.w);
}

__forceinline__ __device__ void LR_atomicAddXYZ(c_number4 *dst, const c_number4 &delta) {
	atomicAdd(&(dst->x), delta.x);
	atomicAdd(&(dst->y), delta.y);
	atomicAdd(&(dst->z), delta.z);
}

__forceinline__ __device__ void LR_atomicAddST(CUDAStressTensor *dst, const CUDAStressTensor &delta) {
	atomicAdd(&(dst->e[0]), delta.e[0]);
	atomicAdd(&(dst->e[1]), delta.e[1]);
	atomicAdd(&(dst->e[2]), delta.e[2]);
	atomicAdd(&(dst->e[3]), delta.e[3]);
	atomicAdd(&(dst->e[4]), delta.e[4]);
	atomicAdd(&(dst->e[5]), delta.e[5]);
}

/**
 * @brief returns the btype (i.e. the fake nucleotide type used to make only certain pairs of particle interacting through HB)
 */
__forceinline__ __device__ int get_particle_btype(const c_number4 &r_i) {
	return __float_as_int(r_i.w) >> 22;
}

__forceinline__ __device__ c_number quad_distance(const c_number4 &r_i, const c_number4 &r_j) {
	const c_number dx = r_j.x - r_i.x;
	const c_number dy = r_j.y - r_i.y;
	const c_number dz = r_j.z - r_i.z;

	return dx * dx + dy * dy + dz * dz;
}

__forceinline__ __device__ int get_particle_type(const c_number4 &r_i) {
	int my_btype = __float_as_int(r_i.w) >> 22;
	return (my_btype > 0) ? (my_btype & 3) : 3 - ((3 - (my_btype)) & 3); // a & 3 is equivalent to a % 4, but much faster
}

__forceinline__ __device__ int get_particle_index(const c_number4 &r_i) {
	int msk = -1 << 22;
	return __float_as_int(r_i.w) & (~msk);
}

__forceinline__ __host__ int get_particle_index_host(const c_number4 &r_i) {
	union {float a; int b;} u;
	u.a = r_i.w;
	int msk = -1 << 22;
	return u.b & (~msk);
}

__forceinline__ __device__ LR_double4 make_LR_double4(const float4 &v) {
	LR_double4 ret;
	ret.x = (double) v.x;
	ret.y = (double) v.y;
	ret.z = (double) v.z;
	ret.w = v.w;

	return ret;
}

__forceinline__ __host__ __device__ c_number4 make_c_number4(const c_number x, const c_number y, const c_number z, const c_number w) {
	c_number4 ret;
	ret.x = x;
	ret.y = y;
	ret.z = z;
	ret.w = (float) w;
	return ret;
}

__forceinline__ __device__ c_number _module(const c_number4 &v) {
	return sqrtf(SQR(v.x) + SQR(v.y) + SQR(v.z));
}

__forceinline__ __device__ c_number _module(const float3 &v) {
	return sqrtf(SQR(v.x) + SQR(v.y) + SQR(v.z));
}

// Necessary to for calculating the torque without storing a separate GPU_matrix on the GPU. Since we have the a1, a2, and a3 vectors anyway, I don't think this is costly. This step might be avoidable if torque and angular momentum were also calculated and stored as quaternions.
__forceinline__ __device__ c_number4 _vectors_c_number4_product(const c_number4 a1, const c_number4 a2, const c_number4 a3, const c_number4 v) {
	c_number4 res = { a1.x * v.x + a2.x * v.y + a3.x * v.z, a1.y * v.x + a2.y * v.y + a3.y * v.z, a1.z * v.x + a2.z * v.y + a3.z * v.z, v.w };

	return res;
}

// Necessary to for calculating the torque without storing a separate GPU_matrix on the GPU. Since we have the a1, a2, and a3 vectors anyway, I don't think this is costly. This step might be avoidable if torque and angular momentum were also calculated and stored as quaternions.
__forceinline__ __device__ c_number4 _vectors_transpose_c_number4_product(const c_number4 a1, const c_number4 a2, const c_number4 a3, const c_number4 v) {
	c_number4 res = { a1.x * v.x + a1.y * v.y + a1.z * v.z, a2.x * v.x + a2.y * v.y + a2.z * v.z, a3.x * v.x + a3.y * v.y + a3.z * v.z, v.w };

	return res;
}

__forceinline__ __device__ c_number4 _cross(const c_number4 v, const c_number4 w) {
	return make_c_number4(v.y * w.z - v.z * w.y, v.z * w.x - v.x * w.z, v.x * w.y - v.y * w.x, (c_number) 0);
}

// LR_DOUBLE4
__forceinline__ __host__ __device__ LR_double4 operator*(LR_double4 v, double c) {
	v.x *= c;
	v.y *= c;
	v.z *= c;
	v.w *= c;
	return v;
}

__forceinline__ __host__ __device__ LR_double4 operator/(LR_double4 v, double c) {
	double inv = 1. / c;
	return v * inv;
}

__forceinline__ __host__ __device__ LR_double4 operator*(double c, LR_double4 v) {
	v.x *= c;
	v.y *= c;
	v.z *= c;
	v.w *= c;
	return v;
}

__forceinline__ __host__ __device__ LR_double4 operator-(LR_double4 a) {
	a.x = -a.x;
	a.y = -a.y;
	a.z = -a.z;
	a.w = -a.w;
	return a;
}

__forceinline__ __host__ __device__ LR_double4 operator+(LR_double4 a) {
	return a;
}

__forceinline__ __host__ __device__ LR_double4 operator+(LR_double4 a, LR_double4 b) {
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
	a.w += b.w;
	return a;
}

__forceinline__ __host__ __device__ void operator+=(LR_double4 &a, LR_double4 b) {
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
	a.w += b.w;
}

__forceinline__ __host__ __device__ void operator*=(LR_double4 &a, double b) {
	a.x *= b;
	a.y *= b;
	a.z *= b;
	a.w *= b;
}

// these two functions add the fourth component (because it doesn't make sense
// to substract them as they contain the energy)
__forceinline__ __device__ LR_double4 operator-(LR_double4 a, LR_double4 b) {
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
	a.w += b.w;
	return a;
}

__forceinline__ __device__ void operator-=(LR_double4 &a, LR_double4 b) {
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
	a.w += b.w;
}

// FLOAT3
__forceinline__ __device__ float3 operator+(float3 a) {
	return make_float3(a.x, a.y, a.z);
}

__forceinline__ __device__ float3 operator-(float3 a) {
	return make_float3(-a.x, -a.y, -a.z);
}

__forceinline__ __device__ float3 operator*(c_number a, float3 b) {
	return make_float3(a * b.x, a * b.y, a * b.z);
}

__forceinline__ __device__ float3 operator*(float3 b, c_number a) {
	return make_float3(a * b.x, a * b.y, a * b.z);
}

__forceinline__ __device__ float3 operator+(float3 a, float3 b) {
	return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

__forceinline__ __device__ float3 operator-(float3 a, float3 b) {
	return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

__forceinline__ __device__ void operator+=(float3 &a, float3 b) {
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
}

__forceinline__ __device__ void operator-=(float3 &a, float3 b) {
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
}

// FLOAT4
__forceinline__ __host__ __device__ float4 operator*(float4 v, float c) {
	return make_float4(v.x * c, v.y * c, v.z * c, v.w * c);
}

__forceinline__ __host__ __device__ float4 operator/(float4 v, float c) {
	float inv = 1.f / c;
	return v * inv;
}

__forceinline__ __host__ __device__ float4 operator*(float c, float4 v) {
	return make_float4(v.x * c, v.y * c, v.z * c, v.w * c);
}

__forceinline__ __host__ __device__ float4 operator-(float4 a) {
	return make_float4(-a.x, -a.y, -a.z, -a.w);
}

__forceinline__ __host__ __device__ float4 operator+(float4 a) {
	return make_float4(a.x, a.y, a.z, a.w);
}

__forceinline__ __host__ __device__ float4 operator+(float4 a, float4 b) {
	return make_float4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

__forceinline__ __host__ __device__ void operator+=(float4 &a, float4 b) {
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
	a.w += b.w;
}

__forceinline__ __host__ __device__ void operator*=(float4 &a, float b) {
	a.x *= b;
	a.y *= b;
	a.z *= b;
	a.w *= b;
}

// these two functions add the fourth component (because it doesn't make sense
// to substract them as they contain the energy)
__forceinline__ __device__ float4 operator-(float4 a, float4 b) {
	return make_float4(a.x - b.x, a.y - b.y, a.z - b.z, a.w + b.w);
}

__forceinline__ __device__ void operator-=(float4 &a, float4 b) {
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
	a.w += b.w;
}

__forceinline__ __device__ c_number4 stably_normalised(const c_number4 &v) {
	c_number max = fmaxf(fmaxf(fabsf(v.x), fabsf(v.y)), fabsf(v.z));
	c_number4 res = v / max;
	return res / _module(res);
}

#endif /* CUDA_LR_COMMON */
