#include "anneal.h"

#ifdef HAVE_MPI
   #include <mpi.h>
#endif

#include <iostream>
#include <cstdlib>
#include <time.h>
#include <sstream>
#include "utils.h"

using namespace std;

int main (int argc, char **argv)
{
	srand(time(0 ));
	vector<string> sequences;
	vector<string> simdata;

#ifdef HAVE_MPI
	 MPI_Init(&argc,&argv);
#endif

 if (argc != 7)
 {
   cout << "Usage " << argv[0] << " log_file_prefix " <<  " parameter_init_file.txt " << " [H/ST/HST/CROSS/AVG] " << " iteration_no " << " beta_init " << " sequence_file simulation_data_file ...."   << endl;  
   return 0;
 }

  string logfile_prefix(argv[1]);
  string param_load_file(argv[2]);
  string simtype(argv[3]);
  int total_iters = atoi(argv[4]);
  double beta_init = atof(argv[5]);

/*
  for(int i = 6; i < argc; i+=2)
   {
     sequences.push_back(string(argv[i]));
     if(i + 1 >= argc)
     {
       cout << string("Error loading simulation data file. Please make sure you include them in the command line afgter sequence file") << endl;
       return 0;
     }
     simdata.push_back(string(argv[i+1]));
   }
*/




try {

    load_fitting_options(argv[6],sequences,simdata);





 logfile_prefix += string("_") ;


#ifdef HAVE_MPI
 int numprocs;
 int myid;
 MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
 MPI_Comm_rank(MPI_COMM_WORLD,&myid);
 if(myid == 0)
    cout << "This is MPI simulation" << endl;


 Anneal_Simulation sim(myid,numprocs,sequences,simdata,param_load_file,simtype.c_str(),beta_init,0.1, logfile_prefix.c_str(),  cout );


 sim.MPI_simulated_annealing(total_iters,beta_init);
//sim.MPI_simulated_annealing_with_maxdifference(total_iters,beta_init);
// sim.MPI_simulated_annealing_with_yields(total_iters,beta_init);
// srand(4);
// sim.MPI_simulated_annealing_with_widths(total_iters,beta_init);

// sim.MPI_retherm_annealing(total_iters,beta_init);
//sim.MPI_STUN_annealing(total_iters,beta_init);
// sim.MPI_retherm_annealing(total_iters,beta_init, 20.0, 0.1, 100, 20, 1000);

 if(myid == 0)
 {
	 cout << "Simulation finished, current status is: " << endl;
	 sim.show_all_melts(cout);
	// sim.show_widths(cout);
 }
 MPI_Finalize();


#endif

#ifndef HAVE_MPI

 Anneal_Simulation sim(sequences,simdata,param_load_file,simtype.c_str(),beta_init,0.1, logfile_prefix.c_str(),  cout );
// Anneal_Simulation sim(param_load_file,simtype.c_str(),beta_init,0.1, logfile_prefix.c_str(),  cout );
// sim.load_sequences_only(sequences,simdata);
// sim.load_fraction_of_states(simdata,0,1);
// sim.calculate_melts_one_sequence_at_time(simdata);
 //srand(4);
// sim.simulated_annealing_with_widths(total_iters,beta_init);


 //sim.STUN_annealing(total_iters,beta_init,0.8,0.2,20,100);
// sim.retherm_simulated_annealing(total_iters,beta_init);
 sim.simulated_annealing(total_iters,beta_init,0.8,0.2,30);
// sim.simulated_annealing(total_iters,beta_init,0.3,0.05,30);
 cout << "Simulation finished, current status is: " << endl;
 //sim.show_widths(cout);
  sim.show_all_melts(cout);
 //sim.show_yields();
#endif



 //sim.simulated_annealing_with_maxdifference(total_iters,beta_init);


 

 
 }
 catch (string &s)
 {
   cout << "Error occured " << s << endl;

 #ifdef HAVE_MPI
  MPI_Abort(MPI_COMM_WORLD,1);
 #endif

 }

 return 1;
}
