#ifndef ANNEAL
#define ANNEAL

#include "pars.h"
#include "fits.h"

#ifdef HAVE_MPI
   #include <mpi.h>
#endif

class Anneal_Simulation
{
 public:

 double olddifference;
 double newdifference ;
 double perturb_max ;

 double sim_beta;
 //int total_sequences;

 double bestdifference;
 vector<string> sequences;
 vector<string> compsequences;
 vector<double> correct_melts;
 vector<FIT_struct> states;  

 vector<int> sequences_in_files;
 
 //double *correct_melts;

 double *attempt_melts;
 double *current_melts;
 
 ostream& log;

 //for mpi use
 int id;
 int proc_no;

 int accepted_moves;
 string simtype;

 //FIT_struct fit_srt;
 PAR_struct str_new;
 PAR_struct str_act;
 PAR_struct str_best;


 ofstream simulation_info_out;
 ofstream best_params_out;
 bool best_params_changed;

 bool average_only;
 bool wc_only; //only watson-crick fitting

 //helper functions
 void set_sim_type(string sim_type); //set type of parameters thst will be adjusted, and whtether it is average, h, or sta, or cross stacking


 //simulation parameters
 int load_data_from_files(const char *seqfilename, const char *energyfilename);
 int load_sequences(const char *seqnamefile, double beta_0, double beta_start, double beta_step, int beta_array_length, PAR_struct &original_parameters, double box, string& simtype);

//full constructor
 Anneal_Simulation(vector<string> &sequence_files, vector<string> &data_files, string& parameter_file, const char *simtype, double _sim_beta,  double pmax, const char * logfile_prefix,  ostream& log ); 
 
 //this is splitting files processing between multiple function calls
 Anneal_Simulation(string& parameter_file, const char *simtype, double _sim_beta,  double pmax, const char * logfile_prefix,  ostream& log );
 void load_sequences_only(vector<string> &sequence_files, vector<string> &data_files);

 //helper functions
 int load_sequences_without_states(const char *seqfilename, const char *energyfilename);
 void load_fraction_of_states(vector<string> &data_files,int myid, int total_id);


 int  calculate_attempt_melts(PAR_struct &pars, int myid, int total_id);
 void calculate_attempt_yields(PAR_struct &pars, int myid, int total_id);
 void calculate_attempt_widths(PAR_struct& pars,int myid, int total_id);

 void calculate_melts_one_sequence_at_time(vector<string> &data_files);


 //helper functions
 double sum_differences(void);
 double pow4_differences(void);

 int calculate_attempt_melts(PAR_struct &pars);
 void calculate_attempt_yields(PAR_struct &pars);
 void calculate_attempt_widths(PAR_struct& pars);
 void calculate_correct_widths(void);

 double get_max_difference(void);
 double get_max_difference_kelvin(void);
 double get_avg_difference_kelvin(void);
 double get_avg_difference(void);

 //single processor MC
 double accept_move(void);
 int simulation_MC_step(double sim_beta);
 void save_best_results(int iter);


 int simulation_STUN_MC_step(double sim_beta);
 int simulation_MC_step_average(double sim_beta); //perturbs only average parameters

 //simulated annealing with l_1 Hamiltonian
 void simulated_annealing(int steps, double beta_start, double high_acceptance = 0.6, double low_acceptance = 0.3, int adapt_after = 10, int output_log_after = 100);

 void simulated_annealing_average(int steps, double beta_start, double high_acceptance = 0.6, double low_acceptance = 0.3, int adapt_after = 10, int output_log_after = 100);
 void retherm_simulated_annealing(int steps, double beta_start, double high_temp = 5.0, double low_temp = 0.2, int adapt_after = 80, int adapt_step = 5, int output_log_after = 800);

 //simulated annealing with l_inf Hamiltonian
 int simulation_MC_step_with_maxdifference(double sim_beta);
 int simulation_MC_step_with_yields(double sim_beta);
 int simulation_MC_step_with_widths(double sim_beta);
// int simulation_MC_step_with_avg_ratios(double sim_beta);


 void simulated_annealing_with_maxdifference(int steps, double beta_start, double high_acceptance = 0.6, double low_acceptance = 0.3, int adapt_after = 10, int output_log_after = 100);
 void simulated_annealing_with_yields(int steps, double beta_start, double high_acceptance = 0.6, double low_acceptance = 0.3, int adapt_after = 10, int output_log_after = 100);
 void simulated_annealing_with_widths(int steps, double beta_start, double high_acceptance = 0.6, double low_acceptance = 0.3, int adapt_after = 10, int output_log_after = 100);

 void STUN_annealing(int steps, double beta_start, double high_acceptance, double low_acceptance , int adapt_after, int output_log_after);

 
 ostream &show(ostream& out);
 void show_all_melts(ostream& out, double cutoff = 0);
 void show_yields(void);
 void show_widths(ostream& out);




 //MPI part
#ifdef HAVE_MPI
 Anneal_Simulation(int rank, int total_proc,vector<string> &sequence_files, vector<string> &data_files,string& parameter_file, const char *simtype, double _sim_beta,  double pmax, const char * logfile_prefix,  ostream& log );

 //double accept_move(void);
 int MPI_simulation_MC_step(double sim_beta);
 int MPI_simulation_MC_step_maxdifference(double sim_beta);
 int MPI_simulation_STUN_MC_step(double beta);
 int MPI_simulation_MC_step_with_yields(double sim_beta);
 int MPI_simulation_MC_step_with_widths(double sim_beta);
 //int MPI_simulation_MC_perturb_avg_ratio_step(double sim_beta);

 int MPI_receive_melts(void);
 int MPI_send_melts(void);
 int MPI_send_pars_new(void);
 int MPI_receive_pars_new(void);

 void MPI_simulated_annealing(int steps, double beta_start, double high_acceptance = 0.6, double low_acceptance = 0.3, int adapt_after = 10, int output_log_after = 100);


 //start old garbage
 void MPI_simulated_annealing_with_maxdifference(int steps, double beta_start, double high_acceptance = 0.6, double low_acceptance = 0.3, int adapt_after = 10, int output_log_after = 100);
 void MPI_retherm_annealing(int steps, double beta_start, double high_temp = 10.0, double low_temp = 0.1, int adapt_after = 100, int adapt_step = 10, int output_log_after = 1000);
 void MPI_STUN_annealing(int steps, double beta_start, double high_acceptance = 0.6, double low_acceptance = 0.3, int adapt_after = 10, int output_log_after = 100);
 void MPI_simulated_annealing_with_widths(int steps, double beta_start, double high_acceptance = 0.6, double low_acceptance = 0.3, int adapt_after = 10, int output_log_after = 100);
// void MPI_simulated_annealing_with_avg_ratios(int steps, double beta_start, double high_acceptance = 0.6, double low_acceptance = 0.3, int adapt_after = 10, int output_log_after = 100);
 void MPI_simulated_annealing_with_yields(int steps, double beta_start, double high_acceptance = 0.6, double low_acceptance = 0.3, int adapt_after = 10, int output_log_after = 100);
 //end old garbage



#endif

};

 inline ostream &operator << (ostream& out, Anneal_Simulation &asim)
 {
   return asim.show(out);
 }



#endif
