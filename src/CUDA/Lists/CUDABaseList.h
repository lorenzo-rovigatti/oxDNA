/**
 * @file    CUDABaseList.h
 * @date    29/set/2010
 * @author  lorenzo
 *
 *
 */

#ifndef CUDABASELIST_H_
#define CUDABASELIST_H_

#include "../CUDAUtils.h"
#include "../cuda_utils/CUDABox.h"

/**
 * @brief Abstract class for list-based force computing on CUDA
 *
 * Child classes have to implement methods to read the input file, update the lists
 * and clean up at the end of the simulation
 */

class CUDABaseList {
protected:
	bool _use_edge;
	int _N;
	CUDABox *_h_cuda_box, *_d_cuda_box;

public:
	CUDABaseList() :
					_use_edge(false),
					_N(-1),
					_h_cuda_box(nullptr),
					_d_cuda_box(nullptr) {
	}
	;
	virtual ~CUDABaseList() {
	}
	;

	virtual void get_settings(input_file &inp) = 0;
	virtual void init(int N, c_number rcut, CUDABox *h_cuda_box, CUDABox *d_cuda_box) {
		_h_cuda_box = h_cuda_box;
		_d_cuda_box = d_cuda_box;
		_N = N;
	}

	virtual void update(c_number4 *poss, c_number4 *list_poss, LR_bonds *bonds) = 0;

	bool use_edge() {
		return _use_edge;
	}

	virtual void clean() = 0;
};

#endif /* CUDABASELIST_H_ */
