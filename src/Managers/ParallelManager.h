/*
 * ParallelManager.h
 *
 *  Created on: 25 Feb 2015
 *      Author: lorenzo
 */

#ifndef PARALLELMANAGER_H_
#define PARALLELMANAGER_H_

#include "SimManager.h"

/**
 * @brief Manages parallel simulations through MPI.
 */
class ParallelManager: public SimManager {
protected:
	int _mpi_rank;
	int _mpi_size;
public:
	ParallelManager(int argc, char *argv[]);
	virtual ~ParallelManager();

	virtual void load_options();
	virtual void init();
};

#endif /* PARALLELMANAGER_H_ */
