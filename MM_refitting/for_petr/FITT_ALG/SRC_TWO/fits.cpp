#include "fits.h"

#include <cstring>
//-------------------------------------------------------------------
STATE_struct::STATE_struct(double weight, double E0_energy, double beta_0)
{
    for (int i = 0; i < PAR_SIZE; i++)
         int_counts[i] = 0;

    other_energies = E0_energy;
    bonded = 0;

    rescaled_hcst_energy = 0;
    unbiased_state = (1.0 / weight) * exp( E0_energy * beta_0);

}
//-----------------------------------------------------------------
double STATE_struct::get_rescaled_energy(PAR_struct &current_pars,double new_T)
 {

    rescaled_hcst_energy = current_pars.get_scaled_energy(int_counts,new_T);
    //return E0_energy;
    return rescaled_hcst_energy + other_energies ;

 }
//-----------------------------------------------------------------
ostream& STATE_struct::show_state(ostream& out)
 {
  out << "State " << " hst energy_rescaled: " << rescaled_hcst_energy << " and the following interaction contributions: " ;
  for(int i =0; i < PAR_SIZE; i++)
  {
    out << int_counts[i] << " " ;
  } 
  cout << endl;
  return out;
 }

//-------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------
void STATE_struct::fill_hb_contribution(string &seq,string& compseq, int sid1, int index_a, int sid2, int index_b, double energy, PAR_struct& original_pars, double beta_0)
{
 try
 {
  if(energy == 0) //no point in filling in 0 energy
  {
	  return;
  }
  other_energies -= energy;
  string pair = "  ";

  if(sid1 == 0)
    pair[0] = seq[index_a];
  else if(sid1 == 1)
    pair[0] = compseq[index_a];

  if(sid2 == 0)
    pair[1] = seq[index_b];
  else if(sid2 == 1)
    pair[1] = compseq[index_b];

  //HERE we need TO DIVIDE by original strength constant

  int fill_id = PAR_struct::return_corresponding_id(string("HBS"),pair);

 // cout << "Fillling in HB " << energy << " as " << original_pars.get_value(fill_id) << endl;


  if(energy < H_CUTOFF)
	  bonded = 1;

  int_counts[fill_id] += energy / original_pars.get_value(fill_id);
 }
 catch(string &s)
 {
	 cerr << "Error occurred for indices" << index_a << " "  << sid1 << " and " << index_b << " " << sid2 << " and energy " << energy << endl;
	 throw s;
 }

}
//------------------------------------------------------------------
void STATE_struct::fill_st_contribution(string &seq,string& compseq, int sid1, int index_a, int sid2, int index_b, double energy, PAR_struct& original_pars, double beta_0)
{
  other_energies -= energy;
  string pair = "  ";

  if(sid1 == 0)
    pair[0] = seq[index_a];
  else if(sid1 == 1)
    pair[0] = compseq[index_a];

  if(sid2 == 0)
    pair[1] = seq[index_b];
  else if(sid2 == 1)
    pair[1] = compseq[index_b];

  //HERE we need TO DIVIDE by original strength constant

  int fill_id = PAR_struct::return_corresponding_id(string("STS"),pair);

 // cout << "Fillling in ST " << energy << " as " << original_pars.get_value(fill_id) << endl;

  int_counts[fill_id] += energy /(   original_pars.get_value(fill_id) *  (1. + original_pars.get_st_T_dep() / beta_0   )  );

}
//---------------------------------------------------------------------------
//------------------------------------------------------------------
void STATE_struct::fill_cross_contribution(string &seq,string& compseq, int sid1, int index_a, int sid2, int index_b, double energy, PAR_struct& original_pars, double beta_0)
{
  other_energies -= energy;
  string pair = "  ";

  if(sid1 == 0)
    pair[0] = seq[index_a];
  else if(sid1 == 1)
    pair[0] = compseq[index_a];

  if(sid2 == 0)
    pair[1] = seq[index_b];
  else if(sid2 == 1)
    pair[1] = compseq[index_b];

  //HERE we need TO DIVIDE by original strength constant

  int fill_id = PAR_struct::return_corresponding_id(string("CROSS"),pair);

//  cout << "Fillling in CR " << energy << " as " << original_pars.get_value(fill_id) << endl;

  int_counts[fill_id] += energy / original_pars.get_value(fill_id);

}
//------------------------------------------------------------------------------------------------------
int FIT_struct::check_no_of_bonded_states(void)
 {
	  int bonds = 0;

	  for(vector<STATE_struct>::iterator i = states.begin(); i != states.end(); i++  )
	  {
		  if( (*i).is_bonded())
			  bonds++;
	  }

	 return bonds;

 }
//------------------------------------------------------------------------------------------------------
FIT_struct::FIT_struct(string& simtype, string& sequence, string& compsequence_, double _beta_0, double beta_start, double beta_step,PAR_struct& old_pars, int _beta_array_length, double box) : beta_array_length(_beta_array_length), sequence(sequence), beta_0(_beta_0), gen_pars(old_pars)
 {
    compsequence = compsequence_;
    
    /*
    for(int i = 0; i < sequence.length(); i++)
    {
      if(sequence[i] == 'A')
         compsequence[sequence.length()-1-i] = 'T';
      if(sequence[i] == 'T')
         compsequence[sequence.length()-1-i] = 'A';
      if(sequence[i] == 'C')
         compsequence[sequence.length()-1-i] = 'G';
      if(sequence[i] == 'G')
         compsequence[sequence.length()-1-i] = 'C';
    }
    */

    sim_type = simtype;

    boxsize = box;
    recorded_states = 0;
    bonded_states = new double[beta_array_length];
    unbonded_states = new double[beta_array_length];
    f_infs = new double[beta_array_length];
    betas = new double[beta_array_length];
    predicted_yields = new double[beta_array_length];
    extra_betas   = new double[beta_array_length];
    extra_f_infs = new double[beta_array_length];
    T_0 = 1. / beta_0;

    symmetric = 0;
    if(sequence == compsequence)
         symmetric = 1;
   

    for(int i = 0; i < beta_array_length; i++)
    {
      betas[i] = beta_start + (double)i * beta_step;
    }


  //  double molarconcentration = 2.6868994 / (boxsize*boxsize*boxsize) ;
  //  get_correct_beta_melt();
  //  fill_in_predicted_yields(molarconcentration);
 }

//-------------------------------------------------------
void FIT_struct::copy_from(const FIT_struct &b)
{
   delete [] bonded_states;
   delete [] unbonded_states;
   delete [] extra_betas;
   delete [] extra_f_infs;

   delete [] f_infs;
   delete [] betas;
   beta_array_length = b.beta_array_length;
   recorded_states = b.recorded_states;
   sequence = b.sequence;
   compsequence = b.compsequence;
   states = b.states;
   gen_pars = b.gen_pars;
   beta_0 = b.beta_0;
   boxsize = b.boxsize; 
   sim_type = b.sim_type;

   bonded_states = new double[beta_array_length];
   unbonded_states = new double[beta_array_length];
   f_infs = new double[beta_array_length];
   betas = new double[beta_array_length];
   predicted_yields = new double[beta_array_length];
   extra_betas   = new double[beta_array_length];
   extra_f_infs = new double[beta_array_length];
   T_0 = 1. / beta_0;

   for(int i = 0; i < beta_array_length; i++)
   {
	   betas[i] = b.betas[i];
	   f_infs[i] = b.f_infs[i];
	   bonded_states[i] = b.bonded_states[i];
	   unbonded_states[i] = b.unbonded_states[i];
	   predicted_yields[i] = b.predicted_yields[i];
	   extra_f_infs[i] = b.extra_f_infs[i];
	   extra_betas[i] =  b.extra_betas[i];

   }
   
}

//----------------------------------------------------
void FIT_struct::reset_values(void)
 {
   recorded_states = 0;
   for(int i = 0; i < beta_array_length; i++)
   {
     bonded_states[i] = unbonded_states[i] = 0;
   }
 }

//----------------------------------------------------

double FIT_struct::get_beta_melt(void)
 {
  fill_in_f_infs();
  int length= beta_array_length;
 // cout << "f_infs" << endl;

  //double *extra_betas = new double[length];
  //double *extra_f_infs = new double[length];


  extra_betas[0] = betas[0];
  extra_f_infs[0] = f_infs[0];
 

  int fill_index = 0; 
  for(int i = 1; i < length; i++)
  {
    if(f_infs[i] - DBL_EPSILON > extra_f_infs[fill_index])
    {   
     fill_index++;
     extra_betas[fill_index] = betas[i];
     extra_f_infs[fill_index] = f_infs[i];
    }
  }
 /*
 if(1.0 - DBL_EPSILON > extra_f_infs[fill_index])
 {
    extra_betas[fill_index+1] = 11.5;
    extra_f_infs[fill_index+1] = 1;
    fill_index++;
 }
*/
/*
for(int i = 0; i < beta_array_length; i++)
 {

   cout << i << " "  << betas[i] << ":" << f_infs[i] << " \n" ;
 }
 cout << endl;

*/

  /*
  cout << "Interpolating " << sequence << endl;

for(int i = 0; i < fill_index+1; i++)
 {

   cout << i << " "  << extra_betas[i] << " " << extra_f_infs[i] << " \n" ;
 }
 cout << endl;
*/

 int check = 0;
 for(int i =1; i < fill_index+1; i++)
 {
   if(extra_f_infs[i] <= extra_f_infs[i-1])
      check = 1;
 }
 if(extra_f_infs[0] > 0.5)
 {

	 cout << "Error in interpolation for sequence " << sequence << ", the lowest beta has yield > 0.5; Try changing interpolation interval" << endl;
	 return 0;
 }
 if(extra_f_infs[fill_index] < 0.5)
 {

 	 cout << "Error in interpolation for sequence " << sequence << ", the highest beta has yield < 0.5; Try changing interpolation interval" << endl;
 	 return 0;
 }


 if(check ||  fill_index < 5)
 {
  cout << "SERIOUS ERROR, non monotonic sequence" << endl;
  for(int i = 0; i < fill_index+1; i++)
  {

   cout << extra_betas[i] << ":" << extra_f_infs[i] << " \n" ;
  }
  cout << endl;

 }

 double beta_melt; 

 {
         gsl_interp_accel *acc  = gsl_interp_accel_alloc ();

         //gsl_spline *spline  = gsl_spline_alloc (gsl_interp_cspline, fill_index+1);
         gsl_spline *spline  = gsl_spline_alloc (gsl_interp_linear, fill_index+1);

         gsl_spline_init (spline,extra_f_infs,extra_betas,fill_index+1);
    

         beta_melt = gsl_spline_eval(spline,0.5,acc);

         gsl_spline_free (spline);
         gsl_interp_accel_free (acc);
 }

 bool sensibility = check_fit_sensibility(extra_betas,extra_f_infs,fill_index+1,beta_melt);
 //cout << " I INTERPOLATED " << beta_melt << endl;
 //cout << "Is the fit sensible? " << sensibility << endl;

 if(!sensibility)
 {
	 cout << "Error while fitting, did not obtain a sensible fit" << endl;
	 return 0;
 }

 //delete [] extra_betas;
 //delete [] extra_f_infs;

 return beta_melt;



}
//--------------------------------------------------------------------
bool FIT_struct::check_fit_sensibility(double *betas, double *yields, int array_size,double beta_melt)
{
	bool sensible = true;

	if (beta_melt <= 0)
	{
		sensible = false;
		return false;
	}

	for(int i = 0; i < array_size; i++)
	{
	    if (yields[i] > 0.5 && betas[i] < beta_melt )
		{
			sensible = false;
			return sensible;
		}
		if(yields[i] < 0.5 && betas[i] > beta_melt)
		{
			sensible = false;
			return sensible;
		}
	}

	return sensible;
}
//-----------------------------------------------------------------
double FIT_struct::get_width_of_melting_curve(PAR_struct& par, double min, double max)
{
	calculate_new_states(par);
	fill_in_f_infs();
	int length= beta_array_length;

	double *extra_betas = new double[length];
	double *extra_f_infs = new double[length];


	extra_betas[0] = betas[0];
	extra_f_infs[0] = f_infs[0];


	int fill_index = 0;
	for(int i = 1; i < length; i++)
	{
		if(f_infs[i] - DBL_EPSILON > extra_f_infs[fill_index])
		{
			fill_index++;
			extra_betas[fill_index] = betas[i];
			extra_f_infs[fill_index] = f_infs[i];
		}
	}

	int check = 0;
	for(int i =1; i < fill_index+1; i++)
	{
		if(extra_f_infs[i] <= extra_f_infs[i-1])
			check = 1;
	}

	if(check ||  fill_index < 5)
	{
		cout << "SERIOUS ERROR, non monotonic sequence" << endl;
		for(int i = 0; i < fill_index+1; i++)
		{

			cout << extra_betas[i] << ":" << extra_f_infs[i] << " \n" ;
		}
		cout << endl;

	}

	double beta_max;
	double beta_min;

	{
		gsl_interp_accel *acc  = gsl_interp_accel_alloc ();
		gsl_spline *spline  = gsl_spline_alloc (gsl_interp_cspline, fill_index+1);

		gsl_spline_init (spline,extra_f_infs,extra_betas,fill_index+1);


		beta_max = gsl_spline_eval(spline,max,acc);
		beta_min = gsl_spline_eval(spline,min,acc);

		gsl_spline_free (spline);
		gsl_interp_accel_free (acc);
	}

	//cout << " I INTERPOLATED " << beta_melt << endl;

	double width = 3000.0/beta_min - 3000.0/beta_max;

	delete [] extra_betas;
	delete [] extra_f_infs;

	return width;




}
//-------------------------------------------------------------------
int is_complementary_sequence(string sequence)
{
   

   string complementary(sequence);
   char complement;
   for(int i = sequence.length() - 1; i >= 0; i--)
   {
      switch(sequence[i])
      {
		  case 'A' : complement = 'T';
		             break;
		  case 'T' : complement = 'A';
		             break;
		  case 'C' : complement = 'G';
		             break;
		  case 'G' : complement = 'C';
		             break;
                 	   
		  default:  printf(" incorrect base letter: %d %c \n",sequence[i],sequence[i]);
		            exit(1);
		            break;
       }
     complementary[(sequence.length()-1)-i] = complement;
   }

   int symmetric = 0;
   
   if(sequence == complementary)
   {
     symmetric = 1; 
    // cout << "Symmetric" << endl;
   }

  
 
  //cout << complementary << endl;
  
  return symmetric;

 
}
//-----------------------------------------------------------------

void FIT_struct::calculate_new_states(PAR_struct &par)
{
   reset_values();
   double Ep;
   //double bonded_weight = 0;
   //double unbonded_weight = 0;
   int bonded = 0;
   double contribution_weight;
   
   for(vector<STATE_struct>::iterator i = states.begin(); i != states.end(); i++)
   {

     //cout << "State Rescaled energy " << (*i).E0_energy << " to " << Ep << " and " << bonded << " "<< endl;
     //cout << " " ;

	  //no T dependence
	 if(par.get_st_T_dep() == 0.)
	 {
		 bonded = (*i).is_bonded();
		 Ep = (*i).get_rescaled_energy(par,1.0/betas[0]);
#ifdef DBG
		 cout << "No T dependence, rescaled energy is " << Ep << endl;
#endif
	 }

     for(int beta_index = 0; beta_index < beta_array_length; beta_index++)
     {
       if(par.get_st_T_dep() != 0)
       {
    		 	bonded = (*i).is_bonded();
    		 	Ep = (*i).get_rescaled_energy(par,1.0/betas[beta_index]);
#ifdef DBG
		 cout << "we have T dependence, rescaled energy is " << Ep << endl;
#endif
       }
       contribution_weight = exp(- betas[beta_index] * Ep) * (*i).get_unbiased_state();
       
     //  contribution_weight = (*i).get_one_over_weight() * exp(- betas[beta_index] * Ep +  (*i).get_E0B0());
    
        /*if(betas[beta_index] == 9.05)
              cout << "weight " << (*i).get_weight() << " with E0: " << (*i).E0_energy << " rescaled to " <<  Ep << " and " << (*i).get_E0B0() << " "  << " and bonds " << (*i).is_bonded() << " " << (*i).currentAT << " "<< (*i).get_unbiased_state()<< " " << " stcontrib: " << (*i).get_rescaled_stacking(gen_pars,par) << " oldst: " << (*i).get_old_stacking() << endl;*/ 
       if(bonded)
       { 
          bonded_states[beta_index] += contribution_weight;
          //cout << bonded_states[beta_index] << " " ;
       }
       else
       { 
          unbonded_states[beta_index] += contribution_weight;
          //cout << unbonded_states[beta_index] << " " ;
       }
     
     }
     // cout << endl;
   }

  /* for(int beta_index = 0; beta_index < beta_array_length; beta_index++)
   {
	   cout << "For " << 8.010648699116 + beta_index * 0.05 << " we have bonded: " << bonded_states[beta_index] << " and unbonded " << unbonded_states[beta_index] << " with ratio " << bonded_states[beta_index] / unbonded_states[beta_index] << endl;
   }*/
    
}
//-----------------------------------------------------------------
double FIT_struct::get_correct_beta_melt(void)
{

  cout << "This function should not have been called!!" << endl;
  throw string("NOT CALLED");

  if(sim_type != "DUPLEXMELT")
  {
    throw string(string("simulation type not implemented") + sim_type);
  } 

  
  string seq = sequence;
  if (sequence == compsequence)
    symmetric = 1;
  else 
    symmetric = 0;

  int length = seq.length();
  double deltaH = 0;
  double deltaS = 0;

  double molarconcentration = 2.686899657  / (boxsize*boxsize*boxsize) ;
  char temp[3] = {'\0','\0','\0'};

  for(int i = 0; i < length - 1; i++)
  {
     temp[0] = seq[i];
     temp[1] = seq[i+1];
     string pair = string(temp);

     if(pair ==  "AA" ) {  deltaH += -6.82 ; deltaS += -19.0; }
     if (pair == "UU") {  deltaH += -6.82; deltaS += -19.0; }

     if (pair == "UA") {  deltaH += -9.38 ; deltaS += -26.7; }
      
     if (pair == "AU") {  deltaH += -7.69; deltaS += -20.5; }
      
     if (pair == "AC") {  deltaH += -10.44; deltaS += -26.9; }
     if (pair == "GU") {  deltaH += -10.44; deltaS += -26.9; }
      
     if (pair == "UG") {  deltaH += -11.4; deltaS += -29.5; }
     if (pair == "CA") {  deltaH += -11.4; deltaS += -29.5; }
        

     if (pair == "UC") {  deltaH += -10.48; deltaS += -27.1; }
     if (pair == "GA") {  deltaH += -10.48; deltaS += -27.1; }
      
     if (pair == "AG") {  deltaH += -12.44; deltaS += -32.5; }
     if (pair == "CU") {  deltaH += -12.44; deltaS += -32.5; }
      
 
     if (pair == "GC") {  deltaH += -10.64; deltaS +=  -26.7; }
      
     if (pair == "CG") {  deltaH += -14.88; deltaS += -36.9; }

     if (pair == "GG") {  deltaH += -13.39; deltaS += -32.7; }
     if (pair == "CC") {  deltaH += -13.39; deltaS += -32.7; }
       
      
  }
    
    //now the init terms
   if(seq[0] == 'A' || seq[0] == 'U')
   {
	deltaH += 3.72;
	deltaS += 10.5;
   }

  if(seq[length-1] == 'A' || seq[length-1] == 'U')
   {
	deltaH += 3.72;
	deltaS += 10.5;
   }

   //initiation 
   deltaH += 3.61;
   deltaS += -1.5;

   deltaH *= 1000.0;

  double divisor = 2.0;
 
   //symetric correcion
   if(symmetric)
   {
     deltaS += -1.4;
     divisor = 0.5;
   }

 // cout << divisor << " " << boxsize << "DH" << deltaH << " ds " << deltaS << " " <<  (molarconcentration/divisor) << endl;
   
  double Tm = deltaH / (deltaS + GAS_CONSTANT * log(molarconcentration/divisor) );

  dH = deltaH;
  dS = deltaS;
  //cout << "Getting beta for " << seq << " and found " <<  3000.0 / Tm  << endl;
  return 3000.0 / Tm; //return Beta_m in simulation units

}
//--------------------------------------------------------------
//-----------------------------------------------------------------
double FIT_struct::get_correct_width( void)
{

  if(sim_type != "DUPLEXMELT")
  {
    throw string(string("simulation type not implemented") + sim_type);
  }

  //HERE WE ASSUME 5' phosphates are not cut out!
 // double salt_correction = 0.368 * + (sequence.length()-1.0) * log(salt_concentration);

  string seq = sequence;
  if (sequence == compsequence)
    symmetric = 1;
  else
    symmetric = 0;

  int length = seq.length();
  double deltaH = 0;
  double deltaS = 0;

  double molarconcentration = 2.686899657 / (boxsize*boxsize*boxsize) ;
  char temp[3] = {'\0','\0','\0'};

  for(int i = 0; i < length - 1; i++)
   {
      temp[0] = seq[i];
      temp[1] = seq[i+1];
      string pair = string(temp);

      if(pair ==  "AA" ) {  deltaH += -6.82 ; deltaS += -19.0; }
      if (pair == "UU") {  deltaH += -6.82; deltaS += -19.0; }

      if (pair == "UA") {  deltaH += -9.38 ; deltaS += -26.7; }

      if (pair == "AU") {  deltaH += -7.69; deltaS += -20.5; }

      if (pair == "AC") {  deltaH += -10.44; deltaS += -26.9; }
      if (pair == "GT") {  deltaH += -10.44; deltaS += -26.9; }

      if (pair == "TG") {  deltaH += -11.4; deltaS += -29.5; }
      if (pair == "CA") {  deltaH += -11.4; deltaS += -29.5; }


      if (pair == "TC") {  deltaH += -10.48; deltaS += -27.1; }
      if (pair == "GA") {  deltaH += -10.48; deltaS += -27.1; }

      if (pair == "AG") {  deltaH += -12.44; deltaS += -32.5; }
      if (pair == "CT") {  deltaH += -12.44; deltaS += -32.5; }


      if (pair == "GC") {  deltaH += -10.64; deltaS +=  -26.7; }

      if (pair == "CG") {  deltaH += -14.88; deltaS += -36.9; }

      if (pair == "GG") {  deltaH += -13.39; deltaS += -32.7; }
      if (pair == "CC") {  deltaH += -13.39; deltaS += -32.7; }


   }

     //now the init terms
    if(seq[0] == 'A' || seq[0] == 'U')
    {
 	deltaH += 3.72;
 	deltaS += 10.5;
    }

   if(seq[length-1] == 'A' || seq[length-1] == 'U')
    {
 	deltaH += 3.72;
 	deltaS += 10.5;
    }

    //initiation
    deltaH += 3.61;
    deltaS += -1.5;


   deltaH *= 1000.0;

  double divisor = 2.0;

   //symetric correcion
   if(symmetric)
   {
     deltaS += -1.4;
     divisor = 0.5;
   }


 // cout << divisor << " " << boxsize << "DH" << deltaH << " ds " << deltaS << " " <<  (molarconcentration/divisor) << endl;
  double width;
  if(!symmetric)
  {
   double Tmin = deltaH / (deltaS - GAS_CONSTANT * log(20.0/molarconcentration ));
   double Tmax = deltaH / (deltaS - GAS_CONSTANT * log(5.0/(molarconcentration*16.0 )) );
   width =  Tmax - Tmin;
  }
  else //to check!
  {
	  double Tmin = deltaH / (deltaS - GAS_CONSTANT * log(20.0/(4.0*molarconcentration) ));
	  double Tmax = deltaH / (deltaS - GAS_CONSTANT * log(5.0/(4.0*molarconcentration*16.0 )) );
	  width =  Tmax - Tmin;

  }


  dH = deltaH;
  dS = deltaS;
  //cout << "Getting beta for " << seq << " and found " <<  3000.0 / Tm  << endl;
  return width; //return Beta_m in simulation units

}

//-------------------------------------------------------------------
void FIT_struct::add_state(string& state_command)
{
    int nuc_id1; 
    int nuc_id2;
    int sid1;
    int sid2;
    double energy;
    stringstream str(state_command);
    double  weight;
    double E0_energy;
    string type;

    recorded_states++;

    str >> weight >> E0_energy ;

    STATE_struct loaded_state(weight,E0_energy, beta_0);   
   

    str >> type >> sid1 >> nuc_id1 >> sid2 >> nuc_id2 >> energy;
    if(!str)
    {
      throw string ( string("Error while reading a state, it contains no interaction energy. Please check format of state file, I loaded ") + str.str()  );
    }
    while (!str == 0)
    { 
       //cout << "Loaded " << type << " " << sid1 << " " << nuc_id1 << " " << sid2 << " " << nuc_id2  << " "<< energy  << endl;
       if(type == "H")
       {
         loaded_state.fill_hb_contribution(sequence,compsequence,sid1, nuc_id1, sid2,nuc_id2,energy,gen_pars,beta_0);
       }
       else if (type == "ST")
       {
         loaded_state.fill_st_contribution(sequence,compsequence,sid1, nuc_id1, sid2, nuc_id2,energy,gen_pars,beta_0);
       }
       else if (type == "CR")
       {
    	 loaded_state.fill_cross_contribution(sequence,compsequence,sid1, nuc_id1, sid2, nuc_id2,energy,gen_pars,beta_0);
       }
       else
       {
          throw string( string("Error reading the file with different fits, unknown interaction ") + type );
       } 
     
       str >> type >> sid1 >> nuc_id1 >> sid2 >> nuc_id2 >> energy;
   
    }    

    //cout << "State loaded " << endl; 
    states.push_back(loaded_state);
    
    
}
//-------------------------------------------------------------------
void FIT_struct::fill_in_f_infs(void )
{
  double ratio = 0;

  //cout << " filling f_infs " << endl;
  for(int beta_ref =0; beta_ref<  beta_array_length; beta_ref++)
  {
    
  

   if(unbonded_states[beta_ref] == 0)
   {
     throw string("Error, division by zero free states");
   }
   
   ratio =  bonded_states[beta_ref]/ unbonded_states[beta_ref];


   //if(symmetric)
   //  ratio *= 2;


   if(sim_type == "DUPLEXMELT")
   {
	   f_infs[beta_ref] =  1.0 + 1.0/(2.0*ratio)  - sqrt( SQR(1.0 + 1.0/(2.0*ratio)  )    - 1.0) ; //ratio;
	  // cout << "DP " << 8.010648699116 + beta_ref * 0.05 << f_infs[beta_ref] << endl;
   }
   else if(sim_type == "HOMOMELT")
   {
	   f_infs[beta_ref] =  1.0 + 1.0/(4.0*ratio)  - sqrt( SQR(1.0 + 1.0/(4.0*ratio)  )    - 1.0) ; //ratio;
	  // cout << "HM " << 8.010648699116 + beta_ref * 0.05 << f_infs[beta_ref] << endl;
   }
   else if(sim_type == "HAIRPINMELT")
   {
	   f_infs[beta_ref] = ratio /(1.0 + ratio) ; //ratio;

   }
   else
   {
	 throw string(string("Unknown type ") + sim_type);
   }
    //cout << betas[beta_ref] << " " << f_infs[beta_ref] << " " << bonded_states[beta_ref] << " " << unbonded_states[beta_ref] << endl;
  }
}
//------------------------------------------------------------------
ostream& FIT_struct::show(ostream& out)
  {
   for(int i = 0; i < states.size(); i++)
      out << states[i] << endl;
   
  
   return out;
  }
//--------------------------------------------------------------------------------
double FIT_struct::calculate_beta_melt_for_parameters(PAR_struct &par)
  {
   
    calculate_new_states(par);
    //fill_in_f_infs();
    return get_beta_melt();

  }

//-------------------------------------------------------------------------------
void FIT_struct::fill_in_predicted_yields(double C0)
{
    //double C0 = 2.6868994 / (boxsize*boxsize*boxsize);
	get_correct_beta_melt();
	for(int i =0; i < beta_array_length; i++)
	{
	   double T = 3000.0 / betas[i];
	   double F_G = exp(- 1.0/(GAS_CONSTANT * T) * (dH - T*dS) );
	   predicted_yields[i] = (1. + 2.* F_G * C0 - sqrt(1.0 + 4.*F_G*C0)  ) / (2.0*F_G*C0 );
       //    cout << "Calculated yield for beta " << betas[i] << " with " <<  F_G << "and " << dH <<  " " << dS << "as " << predicted_yields[i] << endl;
	}
}
//--------------------------------------------------------------------------------
double FIT_struct::get_distance_from_melting_curve(PAR_struct& par)
{
	 calculate_new_states(par);
	 fill_in_f_infs();
	 double dist = 0;
	 for(int i =0 ; i < beta_array_length; i++)
	 {
		 dist += SQR(predicted_yields[i] - f_infs[i] );
	 }

	 return dist / double(beta_array_length);

}
//--------------------------------------------------------------------------------
void FIT_struct::show_yields_for_parameters(PAR_struct& par)
{
	calculate_new_states(par);
	fill_in_f_infs();
	cout << "#beta yield predicted_yield" << endl;
	for(int i =0 ; i < beta_array_length; i++)
	{
		cout << betas[i] << " " << f_infs[i] << " " << predicted_yields[i] << endl;
	}


}
//---------------------------------------------------------------------------------
double FIT_struct::get_avg_beta_melt(void)
{
 int length = sequence.length();
 double deltaH = -10.783125 * (length-1.0);  //-8.2375 * (length-1.0);
 double deltaS =  -27.8875 * (length-1.0); //-22.0188 * (length-1.0);

 deltaH += 3.61 + 3.72;
 deltaS +=  -1.5 + 10.5;

 deltaH *= 1000.0;



 double molarconcentration =  2.686899657 / (boxsize*boxsize*boxsize) ;
 double divisor = 2.0;

    if(sequence == compsequence)
    {
      divisor = 0.5;
    }

 double Tm = deltaH / (deltaS + GAS_CONSTANT * log(molarconcentration/divisor) );

 return 3000.0 / Tm; //return Beta_m in simulation units
}
//---------------------------------------------------------------------------------
double GetBetaForAverage(const char *seq, double boxsize, int symmetric) //this is wrong
{
  int length = strlen(seq) - 1;
  
  double deltaH = -10.783125 * (length);  //-8.2375 * (length-1.0);
  double deltaS =  -27.8875 * (length); //-22.0188 * (length-1.0);

  deltaH += 3.61 + 3.72;
  deltaS +=  -1.5 + 10.5;

  deltaH *= 1000.0;  

  double molarconcentration =  2.686899657 / (boxsize*boxsize*boxsize) ;
  double divisor = 2.0;
 
   if(symmetric)
   {
     divisor = 4.0;
   }

  double Tm = deltaH / (deltaS + GAS_CONSTANT * log(molarconcentration/divisor) );

  return 3000.0 / Tm; //return Beta_m in simulation units
}

//---------------------------------------------------------------------------------
