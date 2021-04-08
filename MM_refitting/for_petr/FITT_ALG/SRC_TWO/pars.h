
#ifndef PARSH
#define PARSH 

#include <iostream>
#include <cstdlib>

#include <cfloat>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>
#include <fstream>

#define GAS_CONSTANT 1.9858775 

#define SQR(x) ((x)*(x))


//this is by what we multiply delta for cross stacking
#define CROSS_DELTA_P 100.

//indicse indicating where start and stop the parameters
#define HB_START 0
#define HB_END 2
#define ST_START 3
#define ST_END 18
#define CROSS_START 19
#define CROSS_END 28
#define ST_T_PART 29
/*
#define HBS_AT_ID 0
#define HBS_CG_ID 1
#define STS_GC_ID 2
#define STS_CG_ID 3
#define STS_GG_ID 4
#define STS_GA_ID 5
#define STS_AG_ID 6
#define STS_TG_ID 7
#define STS_GT_ID 8
#define STS_AT_ID 9
#define STS_TA_ID 10
#define STS_AA_ID 11
*/


#define PAR_SIZE 30

//the number of hydrogen bonding pairs
//#define MAX_HB 2

using namespace std;




inline double rand_01(void)
{
 return double(rand()) / RAND_MAX; 
}


inline int rand_0N(int N)
{
 return rand() % (N+1);
}


inline double beta_to_kelvin(double beta)
{
  return 3000.0 / beta;
}

struct PAR_struct 
{
  static string parameter_names[PAR_SIZE];
  static int wc_correspondence[ST_END+1];
  double vals[PAR_SIZE];
  int allowed_changes[PAR_SIZE];
   

  static int return_corresponding_id(string interaction_type, string pair);


  PAR_struct(void);


  int check_reasonable_limits(void);

  double get_value(int id);
  double get_st_T_dep(void);


  void perturb_parameter(double delta_max);
  void perturb_parameter(double delta_max,bool avg_pars,bool wc_only = false);
  void perturb_wc_only_parameter(double delta_max);


  void perturb_H_only(double delta_max); //perturbs both AT and CG by the same random amount
  void perturb_ST_only(double delta_max); //pertrubs all stacking interactions (10) by the same random amount
  void perturb_CROSS_only(double delta_max);
  void perturb_T_dep_only(double delta_max);

  ostream &show_parameters(ostream& out);
  istream &load_parameters(istream &in);
  void load_from_file(const char* fname) ;

  /*
  void set_relevant_parameters(string sequence);
  void set_only_H_parameters(string sequence);
  */

  
  void disallow_all_changes(void);

  void allow_H(void);
  void allow_ST(void);
  void allow_CROSS(void);
  void allow_ST_dep(void);

  void allow_H_only(void);
  void allow_ST_only(bool T_dependent=false);
  
  double get_scaled_energy(double* different_energies,double new_T);





  istream& readline(istream& in);
  ostream& printline(ostream& out);
  

  //needs to be deleted
  void perturb_parameter_average(double delta_max);
  void perturb_avg_ratio(double delta_max); //perturbs H_bonding / stacking ratio
  void perturb_ATCG_ratio(double delta_max);

};

inline istream &operator >> (istream& in, PAR_struct &par)
{
  return par.readline(in);
}

inline ostream& operator << (ostream& out, PAR_struct &par)
{
  return par.printline(out);
}


/*
inline double get_anneal_beta(int iteration,int iteration_max)
{
  return 1.0;
}
*/
#endif
