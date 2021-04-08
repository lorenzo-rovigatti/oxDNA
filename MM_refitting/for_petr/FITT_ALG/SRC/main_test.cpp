#include "anneal.h"

#ifdef HAVE_MPI
   #include <mpi.h>
#endif

#include <iostream>
#include <cstdlib>
#include <time.h>
#include <sstream>

using namespace std;

int main (int argc, char **argv)
{
	srand(time(0 ));
//	srand(4);

	vector<string> sequences;
	vector<string> simdata;

#ifdef HAVE_MPI
	 MPI_Init(&argc,&argv);
//#endif

	 if (argc < 2)
	 {
	   cout << "Usage " << argv[0] << "\" log_file_prefix " <<  " parameter_init_file.txt " << " [H/ST/HST] " << " iteration_no " << " beta_init " << " sequence_file simulation_data_file .... \" "   << endl;
	   MPI_Finalize();
	   return 0;
	 }

	 stringstream params;
	 params << argv[1];

	 string logfile_prefix;
	 string param_load_file;
	 string simtype;
	 int total_iters ;
	 double beta_init ;

	 params >> logfile_prefix >> param_load_file >> simtype >> total_iters >> beta_init;

	 string fnameSeq, fnameState;
	 params >> fnameSeq >> fnameState;
	 while( !params == 0 )
	 {
		 cout << fnameSeq << " " << fnameState << endl;
		 sequences.push_back(fnameSeq);
		 simdata.push_back(fnameState);
		 //throw("Error while processing parameters");
		 params >> fnameSeq >> fnameState;
		// cout << "Loaded " << fnameSeq << " " << fnameState << " and " << params.good() << endl ;
	 }


#endif

#ifndef HAVE_MPI

 if (argc < 8)
 {
   cout << "Usage " << argv[0] << " log_file_prefix " <<  " parameter_init_file.txt " << " [H/ST/HST] " << " iteration_no " << " beta_init " << " sequence_file simulation_data_file ...."   << endl;  
   return 0;
 }

  string logfile_prefix(argv[1]);
  string param_load_file(argv[2]);
  string simtype(argv[3]);
  int total_iters = atoi(argv[4]);
  double beta_init = atof(argv[5]);

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


#endif



try {






 logfile_prefix += string("_TEST") ;


#ifdef HAVE_MPI
 int numprocs;
 int myid;
 MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
 MPI_Comm_rank(MPI_COMM_WORLD,&myid);
 if(myid == 0)
    cout << "This is MPI simulation" << endl;
 Anneal_Simulation sim(myid,numprocs,sequences,simdata,param_load_file,simtype.c_str(),beta_init,0.1, logfile_prefix.c_str(),  cout );
 sim.MPI_simulated_annealing(0,beta_init);
 if(myid == 0)
 {
     ofstream results(logfile_prefix.c_str());
     results << "#Loaded parameters from " << param_load_file << endl;
     for(int stav = 0; stav < sequences.size(); stav++)
     {
       results << "#Loaded sequences from " << sequences[stav] << " and states from " << simdata[stav] << endl;
     }
	 cout << "Simulation finished, current status is: " << endl;
	 sim.show_all_melts(results);
         results.close();

 }
 MPI_Finalize();


#endif

#ifndef HAVE_MPI
 Anneal_Simulation sim(sequences,simdata,param_load_file,simtype.c_str(),beta_init,0.02, logfile_prefix.c_str(),  cout );
 //Anneal_Simulation sim(sequences,simdata,param_load_file,simtype.c_str(),beta_init,0.1, logfile_prefix.c_str(),  cout );
// Anneal_Simulation sim(param_load_file,simtype.c_str(),beta_init,0.1, logfile_prefix.c_str(),  cout );
// sim.load_sequences_only(sequences,simdata);
// sim.calculate_melts_one_sequence_at_time(simdata);
 sim.simulated_annealing_average(0,beta_init);
 // sim.simulated_annealing(0,beta_init);



 cout << "#Simulation finished, current status is: " << endl;
 sim.show_all_melts(cout);


 ofstream results(logfile_prefix.c_str());
     results << "#Loaded parameters from " << param_load_file << endl;
     for(int stav = 0; stav < sequences.size(); stav++)
     {
       results << "#Loaded sequences from " << sequences[stav] << " and states from " << simdata[stav] << endl;
     }
 for(int stav = 0; stav < sequences.size(); stav++)
 {
	 cout << "#Loaded sequences from " << sequences[stav] << " and states from " << simdata[stav] << endl;
 }
 // cout << "Simulation finished, current status is: " << endl;
 sim.show_all_melts(results);

 //sim.show_yields();
 results.close();
#endif



 //sim.simulated_annealing_with_maxdifference(total_iters,beta_init);


 

 
 }
 catch (string &s)
 {
   cout << "Error occured " << s << endl;

 #ifdef HAVE_MPI
  MPI_Finalize();
 #endif

 }

 return 1;
}
