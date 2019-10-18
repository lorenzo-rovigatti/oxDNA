/*
 * CUDAListFactory.h
 *
 *  Created on: 18/feb/2013
 *      Author: lorenzo
 */

#ifndef CUDALISTFACTORY_H_
#define CUDALISTFACTORY_H_

#include "CUDABaseList.h"

/**
 * @brief Static factory class. Its only public method builds a {@link CUDABaseList CUDA list}.
 *
 * @verbatim
[CUDA_list = no|verlet (Neighbour lists for CUDA simulations. Defaults to 'no'.)]
@endverbatim
 */
class CUDAListFactory {
private:
public:
	CUDAListFactory() = delete;
	virtual ~CUDAListFactory() = delete;

	static CUDABaseList *make_list(input_file &inp);
};

#endif /* CUDALISTFACTORY_H_ */
