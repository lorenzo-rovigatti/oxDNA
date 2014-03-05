/*
 * CUDASharedVerletList.h
 *
 *  Created on: 26/feb/2014
 *      Author: lorenzo
 */

#ifndef CUDASHAREDVERLETLIST_H_
#define CUDASHAREDVERLETLIST_H_

#include "CUDASimpleVerletList.h"

template<typename number, typename number4>
class CUDASharedVerletList: public CUDASimpleVerletList<number, number4> {
protected:
	CUDA_kernel_cfg _kernel_cfg;

	void _set_kernel_cfg();
	void _init_CUDA_verlet_symbols();
public:
	CUDASharedVerletList();
	virtual ~CUDASharedVerletList();

	void init(int N, number box_side, number rcut);
	void update(number4 *poss, number4 *list_poss, LR_bonds *bonds);
	void clean();
};

#endif /* CUDASHAREDVERLETLIST_H_ */
