/**
 * @file    oxDNA_mpi.cpp
 * @date    Created on: 02/set/2010
 * @author  lorenzo (con grosse mani da parte di tutti)
 *
 * @brief Main file for parallel oxDNA simulation
 */

/**
 * @mainpage oxDNA Main page
 * @section Brief
 *
 * oxDNA is a simulation code that was initially conceived as an implementation of the coarse-grained DNA
 * model introduced by T. E. Ouldridge, J. P. K. Doye and A. A. Louis (http://dx.doi.org/10.1063/1.3552946).
 * It has been since reworked and it is now an extensible simulation+analysis framework.
 *
 * See high-level documentation on web: http://dna.physics.ox.ac.uk
 */

#include <mpi.h>

#include "defs.h"
#include "Managers/ParallelManager.h"

using namespace std;

void print_version() {
	fprintf(stdout, "oxDNA %d.%d.%d by Lorenzo Rovigatti, Flavio Romano, Petr Sulc and Benedict Snodin (c) 2013\n", VERSION_MAJOR, VERSION_MINOR, VERSION_STAGE);
	exit(-1);
}

int main(int argc, char *argv[]) {
	Logger::init();

	MPI_Init(&argc,&argv);
	//cout << "This is an MPI simulation" << endl;
	int myid;
	int proc_no;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_no);
	if(myid != 0) {
		Logger::instance()->disable_log();
		//freopen("/dev/null", "w", stderr);
	}

	try {
		TimingManager::init();
		if(argc < 2) throw oxDNAException("Usage is '%s input_file'", argv[0]);
		if(!strcmp(argv[1], "-v")) print_version();

		ParallelManager mysim(argc, argv);
		mysim.load_options();
		MPI_Barrier(MPI_COMM_WORLD);

		OX_DEBUG("Initializing");
		mysim.init();
		MPI_Barrier(MPI_COMM_WORLD);

		OX_LOG(Logger::LOG_INFO, "SVN CODE VERSION: %s", SVN_VERSION);
		OX_LOG(Logger::LOG_INFO, "COMPILED ON: %s", BUILD_TIME);

		OX_DEBUG("Running");
		mysim.run();

		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Finalize();

		OX_LOG(Logger::LOG_INFO, "END OF THE SIMULATION, everything went OK!");
	}
	catch (oxDNAException &e) {
		OX_LOG(Logger::LOG_ERROR, "%s", e.error());
		return 1;
	}

	TimingManager::clear();
	Logger::clear();

	return 0;
}

