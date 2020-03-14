/*
 * CUDAInteractionFactory.h
 *
 *  Created on: 22/feb/2013
 *      Author: lorenzo
 */

#ifndef CUDAINTERACTIONFACTORY_H_
#define CUDAINTERACTIONFACTORY_H_

#include "CUDABaseInteraction.h"

/**
 * @brief Static factory class. Its only public method builds a {@link CUDABaseInteraction CUDA interaction} as specified in the input file.
 */
class CUDAInteractionFactory {
private:
	CUDAInteractionFactory();
public:
	virtual ~CUDAInteractionFactory();

	/**
	 * @brief Builds the interaction instance.
	 *
	 * @param inp
	 * @return a pointer to the newly created interaction
	 */

	static std::shared_ptr<CUDABaseInteraction> make_interaction(input_file &inp);
};

#endif /* CUDAINTERACTIONFACTORY_H_ */
