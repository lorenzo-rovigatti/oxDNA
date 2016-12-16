/**
 * @file    oxDNA.cpp
 * @date    Created on: 02/set/2010
 * @author  lorenzo (con grosse mani da parte di tutti)
 *
 * @brief Main file for oxDNA simulation
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

#include "defs.h"
#include "Managers/SimManager.h"
#include "Utilities/SignalManager.h"

using namespace std;

void print_version() {
	fprintf(stdout, "oxDNA %d.%d.%d by Lorenzo Rovigatti, Flavio Romano, Petr Sulc and Benedict Snodin (c) 2013\n", VERSION_MAJOR, VERSION_MINOR, VERSION_STAGE);
	exit(-1);
}

int main(int argc, char *argv[]) {
	SimManager *mysim;

	try {
		Logger::init();
		SignalManager::manage_segfault();
		TimingManager::init();

		if(argc < 2) throw oxDNAException("Usage is '%s input_file'", argv[0]);
		if(!strcmp(argv[1], "-v")) print_version();

		mysim = new SimManager(argc, argv);
		mysim->load_options();

		OX_DEBUG("Initializing");
		mysim->init();

		OX_LOG(Logger::LOG_INFO, "SVN CODE VERSION: %s", SVN_VERSION);
		OX_LOG(Logger::LOG_INFO, "COMPILED ON: %s", BUILD_TIME);

		OX_DEBUG("Running");
		mysim->run();

		OX_LOG(Logger::LOG_INFO, "END OF THE SIMULATION, everything went OK!");
	}
	catch (oxDNAException &e) {
		OX_LOG(Logger::LOG_ERROR, "%s", e.error());
		return 1;
	}

	delete mysim;
	TimingManager::clear();
	Logger::clear();

	return 0;
}

