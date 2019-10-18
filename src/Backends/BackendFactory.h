/*
 * BackendFactory.h
 *
 *  Created on: Feb 13, 2013
 *      Author: rovigatti
 */

#ifndef BACKENDFACTORY_H_
#define BACKENDFACTORY_H_

#include "SimBackend.h"

/**
 * @brief Static factory class. It exposes a single static method that builds a {@link SimBackend simulation backend}.
 *
 * This class can not be instantiated. It provides a single method that
 * parses the input file and builds the right simulation backend.
 *
 * @verbatim
backend = CPU|CUDA (simulation backend. Defaults to CPU)
[backend_precision = float|double|mixed (Precision at which calculateions are carried out. The mixed precision is only available on CUDA. Defaults to double.)]
[sim_type = MD|MC|VMMC|FFS_MD (Type of the simulation. Supported types are Molecular Dynamics, Monte Carlo, Virtual Move Monte Carlo and Forward Flux Sampling. The first and last ones are also available on CUDA. Defaults to MD.)]
@endverbatim
 */
class BackendFactory {
public:
	BackendFactory() = delete;
	virtual ~BackendFactory() = delete;

	/**
	 * @brief Builds the backend.
	 *
	 * @param inp
	 * @return a pointer to the newly created backend
	 */
	static std::shared_ptr<SimBackend> make_backend(input_file &inp);
};

#endif /* BACKENDFACTORY_H_ */
