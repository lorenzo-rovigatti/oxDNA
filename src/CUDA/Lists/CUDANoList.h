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

class CUDANoList: public CUDABaseList{
public:
	CUDANoList();
	virtual ~CUDANoList();

	void update(c_number4 *poss, c_number4 *list_poss, LR_bonds *bonds) {}
	void clean() {}

	void get_settings(input_file &inp) override;
};

#endif /* CUDANOLIST_H_ */
