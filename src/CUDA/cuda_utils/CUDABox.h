/*
 * CUDABox.h
 *
 *  Created on: 30/mar/2016
 *      Author: lorenzo
 */

#ifndef CUDABOX_H_
#define CUDABOX_H_

#include "../../Boxes/BaseBox.h"

#include <cfloat>
#include <vector_functions.h>
#include "../CUDAUtils.h"

template<typename number, typename number4>
class CUDABox {
protected:
	number _Lx, _Ly, _Lz;
	bool _cubic;

public:
	__host__ __device__ CUDABox() : _cubic(false) {

	}

	__host__ __device__ CUDABox(const CUDABox<number, number4> &b) {
		_cubic = b._cubic;
		_Lx = b._Lx;
		_Ly = b._Ly;
		_Lz = b._Lz;
	}

	~CUDABox() {

	}

	__forceinline__ __device__ number4 minimum_image(number4 &r_i, number4 &r_j) {
		number4 res;
		res.x = r_j.x - r_i.x;
		res.y = r_j.y - r_i.y;
		res.z = r_j.z - r_i.z;

		res.x -= rintf(res.x / _Lx) * _Lx;
		res.y -= rintf(res.y / _Ly) * _Ly;
		res.z -= rintf(res.z / _Lz) * _Lz;

		return res;
	}

	__forceinline__ __device__ number sqr_minimum_image(number4 &r_i, number4 &r_j) {
		number4 mi = minimum_image(r_i, r_j);
		return SQR(mi.x) + SQR(mi.y) + SQR(mi.z);
	}

	__forceinline__ __device__ int compute_cell_index(int N_cells_side[3], float4 &r) {
		int cx = ((r.x/_Lx - floorf(r.x/_Lx)) * (1.f - FLT_EPSILON))*N_cells_side[0];
		int cy = ((r.y/_Ly - floorf(r.y/_Ly)) * (1.f - FLT_EPSILON))*N_cells_side[1];
		int cz = ((r.z/_Lz - floorf(r.z/_Lz)) * (1.f - FLT_EPSILON))*N_cells_side[2];

		return (cz*N_cells_side[1] + cy)*N_cells_side[0] + cx;
	}

	__forceinline__ __device__ int3 compute_cell_spl_idx(int N_cells_side[3], float4 &r) {
		int cx = (r.x/_Lx - floorf(r.x/_Lx)) * (1.f - FLT_EPSILON)*N_cells_side[0];
		int cy = (r.y/_Ly - floorf(r.y/_Ly)) * (1.f - FLT_EPSILON)*N_cells_side[1];
		int cz = (r.z/_Lz - floorf(r.z/_Lz)) * (1.f - FLT_EPSILON)*N_cells_side[2];

		return make_int3(cx, cy, cz);
	}

	__forceinline__ __device__ int compute_cell_index(int N_cells_side[3], LR_double4 &r) {
		int cx = (r.x/_Lx - floor(r.x/_Lx)) * (1. - DBL_EPSILON)*N_cells_side[0];
		int cy = (r.y/_Ly - floor(r.y/_Ly)) * (1. - DBL_EPSILON)*N_cells_side[1];
		int cz = (r.z/_Lz - floor(r.z/_Lz)) * (1. - DBL_EPSILON)*N_cells_side[2];

		return (cz*N_cells_side[1] + cy)*N_cells_side[0] + cx;
	}

	__forceinline__ __device__ int3 compute_cell_spl_idx(int N_cells_side[3], LR_double4 &r) {
		int cx = (r.x/_Lx - floor(r.x/_Lx)) * (1. - DBL_EPSILON)*N_cells_side[0];
		int cy = (r.y/_Ly - floor(r.y/_Ly)) * (1. - DBL_EPSILON)*N_cells_side[1];
		int cz = (r.z/_Lz - floor(r.z/_Lz)) * (1. - DBL_EPSILON)*N_cells_side[2];

		return make_int3(cx, cy, cz);
	}

	/**
	 * @brief (Re-)initialise the object from its CPU counterpart. For now it assumes that the box is orthogonal (cubic or non-cubic)
	 *
	 * @param box
	 */
	void set_CUDA_from_CPU(BaseBox<number> *box) {
		LR_vector<number> sides = box->box_sides();
		change_sides(sides.x, sides.y, sides.z);
	}

	/**
	 * @brief (Re-)initialise the CPU object using the values stored in this object. For now it assumes that the box is orthogonal (cubic or non-cubic)
	 *
	 * @param box
	 */
	void set_CPU_from_CUDA(BaseBox<number> *box) {
		box->init(_Lx, _Ly, _Lz);
	}

	__host__ __device__ number4 box_sides() {
		number4 res;
		res.x = _Lx;
		res.y = _Ly;
		res.z = _Lz;
		res.w = 0.;
		return res;
	}

	void change_sides(number nLx, number nLy, number nLz) {
		_Lx = nLx;
		_Ly = nLy;
		_Lz = nLz;
		if(nLx == nLy && nLy == nLz) _cubic = true;
	}

	number V() {
		return _Lx*_Ly*_Lz;
	}
};

#endif /* CUDABOX_H_ */
