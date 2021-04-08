#ifndef FITSH

#define FITSH

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "pars.h"


#define H_CUTOFF (-0.01f)

double GetBetaForAverage(const char *seq, double boxsize = 20, int symmetric = 0);
//double GetBetaForSequence(const char *seq,double boxsize=20,int symmetric = 0);
int is_complementary_sequence(string sequence);



struct STATE_struct
{
  double int_counts[PAR_SIZE];
  
 // double weight; //weight associated with given state in the simulation
 // double E0_energy; //total_energy_of_the_system
  double other_energies; //total energy of the system - h_Bond_enegy - stacking_energy - cross_stacking_energy
 // double beta_0; //at which beta the simulation was run
  unsigned char bonded;
  double unbiased_state;
  double rescaled_hcst_energy;

 
  STATE_struct(double weight, double E0_energy, double beta_0);

  unsigned char is_bonded(void) { return bonded;}
  double *get_energies(void) { return int_counts;}

  void fill_st_contribution(string &seq,string& compsequence, int sid1, int index_a, int sid2, int index_b, double energy, PAR_struct &original_parameters, double beta_0);
  void fill_hb_contribution(string &seq,string& compsequence, int sid1, int index_a, int sid2, int index_b, double energy, PAR_struct &original_parameters, double beta_0);
  void fill_cross_contribution(string &seq,string& compsequence, int sid1, int index_a, int sid2, int index_b, double energy, PAR_struct &original_parameters, double beta_0);

  double get_unbiased_state(void)  {return unbiased_state;}
  double get_rescaled_energy(PAR_struct &current_pars,double new_T);

  ostream& show_state(ostream& out);


};

inline ostream& operator << (ostream &out, STATE_struct& st)
{
  return st.show_state(out);
}

template <typename T> void really_clear( T & t ) {
    T tmp;
    t.swap( tmp );
}

//FIT class ---------------------------------------------------------------
//-------------------------------------------------------------------------
struct FIT_struct
{
  
 double* bonded_states;
 double* unbonded_states;
 double* f_infs;
 double* betas;

 double *extra_betas; //= new double[length];
 double *extra_f_infs;// = new double[length];

 string sim_type; // = "DUPLEXMELT", "HAIRPINMELT", "HOMOMELT"
 int beta_array_length;
 double *predicted_yields;
 long int recorded_states;

 string sequence;
 string compsequence;
 double dH;
 double dS;

 vector<STATE_struct> states;

 PAR_struct gen_pars;

 double beta_0;
 double T_0; // = 1 / beta_0

 double boxsize; 
 int symmetric;


  FIT_struct(string& simtype, string& sequence, string &compsequence, double beta_0, double beta_start, double beta_step, PAR_struct& old_pars, int beta_array_length , double box=20 );
  
  void reset_values(void); //sets all to 0

  FIT_struct(const FIT_struct &b)
  {
    bonded_states = unbonded_states = f_infs = betas =  extra_f_infs = extra_betas = 0;
    copy_from(b);
  }

  void copy_from(const FIT_struct &b);

  FIT_struct& operator = (const FIT_struct& b)
  {
     copy_from(b);
     return (*this);
  }

 
  int check_no_of_bonded_states(void);

  int get_states_number(void) { return states.size(); }
  void add_state(string& state_command);

  void calculate_new_states(PAR_struct &par);

  void fill_in_f_infs(void);
  double get_beta_melt(void);
  bool check_fit_sensibility(double *betas, double *yields, int array_size,double beta_melt); //checks if fit is reasonable

  double calculate_beta_melt_for_parameters(PAR_struct &par);

  double get_correct_beta_melt(void); //data from Turner (99)
  double get_avg_beta_melt(void); //Turner (99), average base dH and dS

  double get_correct_width(void); //should be 0.5
  void clear(void)  {  really_clear(states); }
  ostream& show(ostream& out);
 
  ~FIT_struct() 
  {delete [] bonded_states; delete [] unbonded_states; delete [] betas; delete [] f_infs; delete [] predicted_yields;  delete [] extra_betas;
   delete [] extra_f_infs;}


  //those need to be rewritten
   double get_width_of_melting_curve(PAR_struct& par, double min=0.2, double max=0.8);
   double get_distance_from_melting_curve(PAR_struct& par);
   void fill_in_predicted_yields(double C0);
   void show_yields_for_parameters(PAR_struct& par);

};



inline ostream& operator << (ostream& out, FIT_struct& fstr)
{
 return fstr.show(out);
}





#endif
