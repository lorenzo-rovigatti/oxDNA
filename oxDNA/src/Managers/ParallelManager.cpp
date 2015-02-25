/*
 * ParallelManager.cpp
 *
 *  Created on: 25 Feb 2015
 *      Author: lorenzo
 */

#include <mpi/mpi.h>

#include "ParallelManager.h"

ParallelManager::ParallelManager(int argc, char *argv[]) : SimManager(argc, argv) {
	MPI_Comm_rank(MPI_COMM_WORLD, &_mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &_mpi_size);
}

ParallelManager::~ParallelManager() {

}

void ParallelManager::load_options() {
	std::string new_prefix = Utils::sformat("output_prefix = mpi_%d_", _mpi_rank);
	addInput(&_input, new_prefix);

	SimManager::load_options();
}

void ParallelManager::init() {
	// a different seed for each process
	_seed = _seed + _mpi_rank;
	SimManager::init();
}
