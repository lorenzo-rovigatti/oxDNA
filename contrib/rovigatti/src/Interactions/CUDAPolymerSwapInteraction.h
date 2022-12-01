/*
 * CUDAPolymerSwapInteraction.h
 *
 *  Created on: 24/mar/2020
 *      Author: lorenzo
 */

#ifndef CUDAPOLYMERSWAPINTERACTION_H_
#define CUDAPOLYMERSWAPINTERACTION_H_

#include "CUDA/Interactions/CUDABaseInteraction.h"
#include "PolymerSwapInteraction.h"

__device__ __host__ __align__(16) struct CUDAStressTensor {
	c_number e[6];

	__device__ __host__ CUDAStressTensor() : e{0} {

	}

	__device__ __host__ CUDAStressTensor(c_number e0, c_number e1, c_number e2, c_number e3, c_number e4, c_number e5) :
		e{e0, e1, e2, e3, e4, e5} {

	}

	__device__ __host__ inline CUDAStressTensor operator+(const CUDAStressTensor &other) const {
		return CUDAStressTensor(
				e[0] + other.e[0],
				e[1] + other.e[1],
				e[2] + other.e[2],
				e[3] + other.e[3],
				e[4] + other.e[4],
				e[5] + other.e[5]);
	}

	__device__ __host__ inline void operator+=(const CUDAStressTensor &other) {
		e[0] += other.e[0];
		e[1] += other.e[1];
		e[2] += other.e[2];
		e[3] += other.e[3];
		e[4] += other.e[4];
		e[5] += other.e[5];
	}

	__host__ inline void operator*=(const c_number other) {
		e[0] *= other;
		e[1] *= other;
		e[2] *= other;
		e[3] *= other;
		e[4] *= other;
		e[5] *= other;
	}

	__host__ inline void operator/=(const c_number other) {
		e[0] /= other;
		e[1] /= other;
		e[2] /= other;
		e[3] /= other;
		e[4] /= other;
		e[5] /= other;
	}
};

/**
 * @brief CUDA implementation of the {@link PolymerSwapInteraction interaction}.
 */
class CUDAPolymerSwapInteraction: public CUDABaseInteraction, public PolymerSwapInteraction {
private:
	c_number4 *_d_three_body_forces = nullptr;
	int *_d_bonded_neighs = nullptr;
	CUDAStressTensor *_d_st = nullptr, *_h_st = nullptr;

public:
	CUDAPolymerSwapInteraction();
	virtual ~CUDAPolymerSwapInteraction();

	void get_settings(input_file &inp);
	void cuda_init(int N);
	c_number get_cuda_rcut() {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box);
};

extern "C" BaseInteraction *make_CUDAPolymerSwapInteraction() {
	return new CUDAPolymerSwapInteraction();
}

#endif /* CUDAPOLYMERSWAPINTERACTION_H_ */
