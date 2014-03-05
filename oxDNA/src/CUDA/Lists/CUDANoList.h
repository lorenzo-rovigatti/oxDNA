/**
 * @file    CUDANoList.h
 * @date    29/set/2010
 * @author  lorenzo
 *
 *
 */

#ifndef CUDANOLIST_H_
#define CUDANOLIST_H_

#include "CUDABaseList.h"

/**
 * @brief Implements a O(N^2) type of simulation (each particle interact with each other) with CUDA.
 */
template<typename number, typename number4>
class CUDANoList: public CUDABaseList<number, number4> {
public:
	CUDANoList();
	virtual ~CUDANoList();

	void init(int N, number box_side, number rcut) {}
	void update(number4 *poss, number4 *list_poss, LR_bonds *bonds) {}
	void clean() {}

	void get_settings(input_file &inp) {}
};

#endif /* CUDANOLIST_H_ */
