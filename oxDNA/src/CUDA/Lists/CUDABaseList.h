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

/**
 * @brief Abstract class for list-based force computing on CUDA
 *
 * Child classes have to implement methods to read the input file, update the lists
 * and clean up at the end of the simulation
 */
template<typename number, typename number4>
class CUDABaseList {
protected:
	int _N;
	number _box_side;
	bool _use_edge;

public:
	CUDABaseList() : _use_edge(false) {};
	virtual ~CUDABaseList() {};

	virtual void get_settings(input_file &inp) = 0;
	virtual void init(int N, number box_side, number rcut) {
		_N = N;
		_box_side = box_side;
	}

	virtual void update(number4 *poss, number4 *list_poss, LR_bonds *bonds) = 0;
	bool use_edge() { return _use_edge; }

	virtual void clean() = 0;
};

#endif /* CUDABASELIST_H_ */
