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

#include <curand_kernel.h>

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

template<typename number4>
__forceinline__ __device__ void LR_atomicAdd(number4 *dst, number4 delta) {
	atomicAdd(&(dst->x), delta.x);
	atomicAdd(&(dst->y), delta.y);
	atomicAdd(&(dst->z), delta.z);
	atomicAdd(&(dst->w), delta.w);
}

template<typename number4>
__forceinline__ __device__ void LR_atomicAddXYZ(number4 *dst, number4 delta) {
	atomicAdd(&(dst->x), delta.x);
	atomicAdd(&(dst->y), delta.y);
	atomicAdd(&(dst->z), delta.z);
}

/**
 * @brief returns the btype (i.e. the fake nucleotide type used to make only certain pairs of particle interacting through HB)
 */
template <typename number, typename number4>
__forceinline__ __device__ int get_particle_btype(const number4 &r_i) {
	return __float_as_int(r_i.w) >> 22;
}

template <typename number, typename number4>
__forceinline__ __device__ number quad_distance(const number4 &r_i, const number4 &r_j) {
	const number dx = r_j.x - r_i.x;
	const number dy = r_j.y - r_i.y;
	const number dz = r_j.z - r_i.z;

	return dx*dx + dy*dy + dz*dz;
}

template <typename number, typename number4>
__forceinline__ __device__ int get_particle_type(const number4 &r_i) {
	int my_btype = __float_as_int(r_i.w) >> 22;
	if (my_btype >= 0 && my_btype <=3) return my_btype;
	if (my_btype > 0) return my_btype % 4;
	else return 3 - ((3 - (my_btype)) % 4);
}

__forceinline__ __device__ LR_double4 make_LR_double4(const float4 &v) {
	LR_double4 ret;
	ret.x = (double) v.x;
	ret.y = (double) v.y;
	ret.z = (double) v.z;
	ret.w = v.w;

	return ret;
}

template <typename number, typename number4>
__forceinline__ __device__ number4 make_number4(const number x, const number y, const number z, const number w) {
	number4 res = {x, y, z, w};
	return res;
}

template <typename number, typename number4>
__forceinline__ __device__ number _module(const number4 v) {
	return sqrtf(SQR(v.x) + SQR(v.y) + SQR(v.z));
}

template <typename number, typename number4>
__forceinline__ __device__ number4 _matrix_number4_product(const LR_GPU_matrix<number> m, const number4 v) {
	number4 res = {
		m.e[0]*v.x + m.e[1]*v.y + m.e[2]*v.z,
		m.e[3]*v.x + m.e[4]*v.y + m.e[5]*v.z,
		m.e[6]*v.x + m.e[7]*v.y + m.e[8]*v.z,
		v.w
	};

	return res;
}

template <typename number, typename number4>
__forceinline__ __device__ number4 _matrix_transpose_number4_product(const LR_GPU_matrix<number> m, const number4 v) {
	number4 res = {
		res.x = m.e[0]*v.x + m.e[3]*v.y + m.e[6]*v.z,
		res.y = m.e[1]*v.x + m.e[4]*v.y + m.e[7]*v.z,
		res.z = m.e[2]*v.x + m.e[5]*v.y + m.e[8]*v.z,
		res.w = v.w
	};

	return res;
}

template <typename number, typename number4>
__forceinline__ __device__ number4 _cross(const number4 v, const number4 w) {
	return make_number4<number, number4>(v.y*w.z - v.z*w.y, v.z*w.x - v.x*w.z, v.x*w.y - v.y*w.x, (number) 0);
}

template <typename number>
__forceinline__ __device__ LR_GPU_matrix<number> _matrix_matrix_product(const LR_GPU_matrix<number> m, const LR_GPU_matrix<number> n) {
	LR_GPU_matrix<number> res = {
		m.e[0]*n.e[0] + m.e[1]*n.e[3] + m.e[2]*n.e[6],
		m.e[0]*n.e[1] + m.e[1]*n.e[4] + m.e[2]*n.e[7],
		m.e[0]*n.e[2] + m.e[1]*n.e[5] + m.e[2]*n.e[8],

		m.e[3]*n.e[0] + m.e[4]*n.e[3] + m.e[5]*n.e[6],
		m.e[3]*n.e[1] + m.e[4]*n.e[4] + m.e[5]*n.e[7],
		m.e[3]*n.e[2] + m.e[4]*n.e[5] + m.e[5]*n.e[8],

		m.e[6]*n.e[0] + m.e[7]*n.e[3] + m.e[8]*n.e[6],
		m.e[6]*n.e[1] + m.e[7]*n.e[4] + m.e[8]*n.e[7],
		m.e[6]*n.e[2] + m.e[7]*n.e[5] + m.e[8]*n.e[8]
	};

	return res;
}

// LR_DOUBLE4

__forceinline__ __device__ LR_double4 operator*(LR_double4 v, double c) {
	return make_number4<double, LR_double4>(v.x*c, v.y*c, v.z*c, v.w);
}

__forceinline__ __device__ LR_double4 operator/(LR_double4 v, double c) {
	double inv = 1./c;
	return v*inv;
}

__forceinline__ __device__ LR_double4 operator*(double c, LR_double4 v) {
	return make_number4<double, LR_double4>(v.x*c, v.y*c, v.z*c, v.w);
}

__forceinline__ __device__ LR_double4 operator-(LR_double4 a) {
	return make_number4<double, LR_double4>(-a.x, -a.y, -a.z, -a.w);
}

__forceinline__ __device__ LR_double4 operator+(LR_double4 a) {
	return make_number4<double, LR_double4>(a.x, a.y, a.z, a.w);
}

__forceinline__ __device__ LR_double4 operator+(LR_double4 a, LR_double4 b) {
	return make_number4<double, LR_double4>(a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w);
}

__forceinline__ __device__ void operator+=(LR_double4 &a, LR_double4 b) {
	a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}

__forceinline__ __device__ void operator*=(LR_double4 &a, double b) {
	a.x *= b; a.y *= b; a.z *= b; a.w *= b;
}

// these two functions add the fourth component (because it doesn't make sense
// to substract them as they contain the energy)
__forceinline__ __device__ LR_double4 operator-(LR_double4 a, LR_double4 b) {
	return make_number4<double, LR_double4>(a.x-b.x, a.y-b.y, a.z-b.z, a.w+b.w);
}

__forceinline__ __device__ void operator-=(LR_double4 &a, LR_double4 b) {
	a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w += b.w;
}


// FLOAT4
__forceinline__ __device__ float4 operator*(float4 v, float c) {
	return make_float4(v.x*c, v.y*c, v.z*c, v.w);
}

__forceinline__ __device__ float4 operator/(float4 v, float c) {
	float inv = 1.f/c;
	return v*inv;
}

__forceinline__ __device__ float4 operator*(float c, float4 v) {
	return make_float4(v.x*c, v.y*c, v.z*c, v.w);
}

__forceinline__ __device__ float4 operator-(float4 a) {
	return make_float4(-a.x, -a.y, -a.z, -a.w);
}

__forceinline__ __device__ float4 operator+(float4 a) {
	return make_float4(a.x, a.y, a.z, a.w);
}

__forceinline__ __device__ float4 operator+(float4 a, float4 b) {
	return make_float4(a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w);
}

__forceinline__ __device__ void operator+=(float4 &a, float4 b) {
	a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}

__forceinline__ __device__ void operator*=(float4 &a, float b) {
	a.x *= b; a.y *= b; a.z *= b; a.w *= b;
}

// these two functions add the fourth component (because it doesn't make sense
// to substract them as they contain the energy)
__forceinline__ __device__ float4 operator-(float4 a, float4 b) {
	return make_float4(a.x-b.x, a.y-b.y, a.z-b.z, a.w+b.w);
}

__forceinline__ __device__ void operator-=(float4 &a, float4 b) {
	a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w += b.w;
}

#endif /* CUDA_LR_COMMON */
