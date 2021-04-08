/*
 * binaryMixture.c
 *
 *  Created on: 02/set/2010
 *      Author: lorenzo (con grosse mani da parte di tutti)
 */

extern "C" {
#include "parse_input/parse_input.h"
}

#include "defs.h"
#include "SimManager.h"

using namespace std;

int main(int argc, char *argv[]) {

	IOManager IO;

#ifdef HAVE_MPI
	MPI_Init(&argc,&argv);
	//cout << "This is an MPI simulation" << endl;
	int myid;
	int proc_no;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&proc_no);
	if(myid != 0) {
		IO.disable_log();
		//freopen("/dev/null", "w", stderr);
	}
#endif

	if(argc < 2) IO.die("Usage is '%s input_file'", argv[0]);
	if(!strcmp(argv[1], "-v")) IO.print_version();

	SimManager mysim(IO, argv[1]);

	IO.debug("Loading options");
	mysim.load_options();
#ifdef HAVE_MPI
	MPI_Barrier (MPI_COMM_WORLD);
#endif

	IO.debug("Initializing");
	mysim.init();
#ifdef HAVE_MPI
	MPI_Barrier (MPI_COMM_WORLD);
#endif

	IO.log(IO.LOG_INFO, "SVN CODE VERSION: %s", SVN_VERSION);
	IO.log(IO.LOG_INFO, "COMPILED ON: %s", BUILD_TIME);
	IO.log(IO.LOG_INFO, "MODEL VERSION: %s", MODEL_VERSION);

	IO.debug("Running");
	mysim.run();

#ifdef HAVE_MPI
	MPI_Barrier (MPI_COMM_WORLD);
#endif

	IO.debug("Cleaning");
	mysim.clean();

#ifdef HAVE_MPI
	MPI_Finalize();
#endif

	return 0;
}

