/**
 * @file    DNAnalysis.cpp
 *
 * @brief Main file for the analysis of oxDNA trajectories
 */

/**
 * DNAnalysis is an off-line tool used to analyse oxDNA trajectories.
 *
 * See high-level documentation on web: http://dna.physics.ox.ac.uk
 */

#include "defs.h"
#include "Managers/AnalysisManager.h"
#include "Utilities/SignalManager.h"
#include "Utilities/oxDNAException.h"
#include "Utilities/Timings.h"
#include "Utilities/parse_input/parse_input.h"

using namespace std;

void print_version() {
	fprintf(stdout, "Configuration analyser for oxDNA %s by Lorenzo Rovigatti, Flavio Romano, Petr Sulc and Benedict Snodin (c) 2013\n", RELEASE);
	exit(-1);
}

int main(int argc, char *argv[]) {
	try {
		Logger::init();
		SignalManager::manage_segfault();
		TimingManager::init();

		if(argc < 2) {
			throw oxDNAException("Usage is '%s input_file'", argv[0]);
		}
		if(!strcmp(argv[1], "-v")) {
			print_version();
		}

		input_file input(true);
		input.init_from_command_line_args(argc, argv);

		AnalysisManager myanalysis(input);
		myanalysis.load_options();

		OX_DEBUG("Initializing");
		myanalysis.init();

		OX_LOG(Logger::LOG_INFO, "RELEASE: %s", RELEASE);
		OX_LOG(Logger::LOG_INFO, "GIT COMMIT: %s", GIT_COMMIT);
		OX_LOG(Logger::LOG_INFO, "COMPILED ON: %s", BUILD_TIME);

		OX_DEBUG("Running");
		OX_LOG(Logger::LOG_INFO, "Starting analysis");
		myanalysis.analysis();

		OX_LOG(Logger::LOG_INFO, "END OF THE ANALYSIS, everything went OK!");
	}
	catch (oxDNAException &e) {
		OX_LOG(Logger::LOG_ERROR, "%s", e.what());
		return 1;
	}

	TimingManager::clear();

	return 0;
}
