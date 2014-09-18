/*
 * delta.cpp
 *
 *  Created on: 20/set/2011
 *      Author: lorenzo
 */

extern "C" {
#include "../../parse_input/parse_input.h"
}

#include "../../defs.h"
#include "DeltaManager.h"

using namespace std;

int main(int argc, char *argv[]) {
	IOManager IO;

	if(argc < 2) IO.die("Usage is '%s input_file'", argv[0]);
	if(!strcmp(argv[1], "-v")) IO.print_version();

	DeltaManager mysim(IO, argv[1]);

	IO.debug("Loading options");
	mysim.load_options();

	IO.debug("Initializing");
	mysim.init();

	IO.log(IO.LOG_INFO, "SVN CODE VERSION: %s", SVN_VERSION);
	IO.log(IO.LOG_INFO, "COMPILED ON: %s", BUILD_TIME);
	IO.log(IO.LOG_INFO, "MODEL VERSION: %s", MODEL_VERSION);

	IO.debug("Running");
	mysim.run();

	IO.debug("Cleaning");
	mysim.clean();

	return 0;
}