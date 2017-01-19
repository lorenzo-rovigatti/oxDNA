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

//This is the most commonly called quaternion to matrix conversion. 
template<typename number, typename number4>
__device__ void get_vectors_from_quat(GPU_quat<number> &q, number4 &a1, number4 &a2, number4 &a3) {
	number sqx = q.x*q.x;
	number sqy = q.y*q.y;
	number sqz = q.z*q.z;
	number sqw = q.w*q.w;
	number xy = q.x*q.y;
	number xz = q.x*q.z;
	number xw = q.x*q.w;
	number yz = q.y*q.z;
	number yw = q.y*q.w;
	number zw = q.z*q.w;

	a1.x = (sqx-sqy-sqz+sqw);
	a2.x = 2*(xy-zw);
	a3.x = 2*(xz+yw);
	a1.y = 2*(xy+zw);
	a2.y = (-sqx+sqy-sqz+sqw);
	a3.y = 2*(yz-xw);
	a1.z = 2*(xz-yw);
	a2.z = 2*(yz+xw);
	a3.z = (-sqx-sqy+sqz+sqw);
}

template <typename number>
__forceinline__ __device__ GPU_quat<number> quat_multiply(GPU_quat<number> &a, GPU_quat<number> &b) {
	GPU_quat<number> p;
	p.w = a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z;
	p.x = a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y;
	p.y = a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x;
	p.z = a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w;
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

template <typename number, typename number4>
__forceinline__ __device__ int get_particle_index(const number4 &r_i) {
	int msk = -1 << 22;
	return __float_as_int(r_i.w) & (~msk);
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


//Necessary to for calculating the torque without storing a seperate GPU_matrix on the GPU. Since we have the a1, a2, and a3 vectors anyway, I don't think this is costly. This step might be avoidable if torque and angular momentum were also calculated and stored as quaternions.
template <typename number4>
__forceinline__ __device__ number4 _vectors_number4_product(const number4 a1, const number4 a2, const number4 a3, const number4 v) {
	number4 res = {
		a1.x*v.x + a2.x*v.y + a3.x*v.z,
		a1.y*v.x + a2.y*v.y + a3.y*v.z,
		a1.z*v.x + a2.z*v.y + a3.z*v.z,
		v.w
	};
	
	return res;
}


//Necessary to for calculating the torque without storing a seperate GPU_matrix on the GPU. Since we have the a1, a2, and a3 vectors anyway, I don't think this is costly. This step might be avoidable if torque and angular momentum were also calculated and stored as quaternions.
template <typename number4>
__forceinline__ __device__ number4 _vectors_transpose_number4_product(const number4 a1, const number4 a2, const number4 a3, const number4 v) {
	number4 res = {
		a1.x*v.x + a1.y*v.y + a1.z*v.z,
		a2.x*v.x + a2.y*v.y + a2.z*v.z,
		a3.x*v.x + a3.y*v.y + a3.z*v.z,
		v.w
	};

	return res;
}

template <typename number, typename number4>
__forceinline__ __device__ number4 _cross(const number4 v, const number4 w) {
	return make_number4<number, number4>(v.y*w.z - v.z*w.y, v.z*w.x - v.x*w.z, v.x*w.y - v.y*w.x, (number) 0);
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
