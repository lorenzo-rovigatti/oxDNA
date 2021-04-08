#include "anneal.h"

//---------------------------------------------------------------------------------
int Anneal_Simulation::load_sequences(const char *seqnamefile, double beta_0, double beta_start, double beta_step, int beta_array_length, PAR_struct &original_parameters, double box, string& simtype)
{
 //sequences.resize(1000);

 
 ifstream fin(seqnamefile);
 int lines = 0;
 char readline[4096];
 string sequence;
 string compsequence;
 double ideal_T; //this is the temperature at which it shoudl be melting, in kelvin

// line >> beta_start >> beta_step >> beta_array_length;
 if(id == 0)
     cout << " Loading sequences from file " << seqnamefile << endl;

 fin.getline(readline,4096);
 while(fin.good() )
 {
  string line(readline);
  if(line.length() > 1 && line[0] != ' ' && line[0] != '\n')
  {
    stringstream str(line);
    str  >> sequence >> compsequence >> ideal_T;

#ifdef DBG
    cout << "Adding " << sequence  << " " << compsequence << " " << ideal_T << endl;
#endif

    FIT_struct fitstr(simtype,sequence,compsequence,beta_0,beta_start,beta_step,original_parameters,beta_array_length, box);
    sequences.push_back(sequence);
    compsequences.push_back(compsequence);
#ifndef WITHWIDTHS
    correct_melts.push_back(3000./ideal_T); //we store the correct beta in simulation units: THIS FOR MELTING
#else
    correct_melts.push_back(ideal_T); //we store the correct beta in simulation units:  THIS IS FOR WIDTH
#endif


    states.push_back(fitstr);
    lines++;
  }
  fin >> line;
 }  

 //total_sequences = lines;
 return lines;

}
//---------------------------------------------------------------------------------
Anneal_Simulation::Anneal_Simulation(vector<string> &sequence_files, vector<string> &data_files, string& parameter_file, const char *simtype, double _sim_beta,  double pmax, const char * logfile_prefix,  ostream& log ) 
 : log(log), simtype(simtype) {
  
   id = 0;
   proc_no = 1;
   string infoname =  string(logfile_prefix)  + string(simtype) + string(".log");
   simulation_info_out.open(infoname.c_str());

   string parname =  string(logfile_prefix)  + string(simtype) +  string(".par");
   best_params_out.open(parname.c_str());
   best_params_changed = false;
  
   sim_beta = _sim_beta;
   olddifference = 100;
   newdifference = 100;
   bestdifference = 100;
   perturb_max = pmax;
   int total_sequences = 0;

   if(sequence_files.size() != data_files.size())
   {
     throw string("Incompatible dimensions ");
   }
  

   log << "There are in total " << sequence_files.size() << endl; 
   for(int i = 0; i < sequence_files.size(); i++)
   {
     log << "Loading sequences  and energy data from " << sequence_files[i] << " and " << data_files[i] << endl;
     total_sequences += load_data_from_files(sequence_files[i].c_str(),data_files[i].c_str());
   } 

   log << "Loaded number of sequences: " << total_sequences << endl;

   //correct_melts = new double[total_sequences];
   attempt_melts = new double[total_sequences];
   current_melts = new double[total_sequences];
  
   if(total_sequences != states.size() )
   {
     throw string("Problem with mismatching number of sequences and fitted strcture. A problem with format of the sequnce file ?");
   }
  

   log << "Loading parameters from " << parameter_file << endl;  
   str_act.load_from_file(parameter_file.c_str());
   log << "loaded " << str_act << endl;
   log << "Setting relevant parameters " << endl;
   for(int i = 0; i < total_sequences; i++)
   { 
	  //CHANGE HERE for SEQUENCE dependent and SEQUENCE-independent
      //correct_melts[i] = states[i].get_correct_beta_melt(); //GetBetaForSequence(sequences[i].c_str(),states[i].get_box_size(),is_complementary_sequence(sequences[i]));

     // correct_melts[i] = GetBetaForAverage(sequences[i].c_str()); //GetBetaForSequence(sequences[i].c_str(),states[i].get_box_size(),is_complementary_sequence(sequences[i]));

     //  str_act.set_relevant_parameters(sequences[i].c_str());
   }

   set_sim_type(simtype);


 // cout << "States bonded:" << endl;
 // states[0].check_no_of_bonded_states();
  
  str_best = str_new = str_act;
   
  log << "Calculating melting temperatures for loaded parameters" << endl;
  for(int i = 0; i < total_sequences; i++)
  {
    current_melts[i] = states[i].calculate_beta_melt_for_parameters(str_act);
    log << sequences[i] << ": " << current_melts[i]  << " " << 3000. / current_melts[i] << endl;
   
  }


 }

//-----------------------------------------------------------------------------------
void Anneal_Simulation::set_sim_type(string simtype)
{
	   average_only = false;

	   bool changed = false;

	   str_act.disallow_all_changes();
	   if( simtype.find(string("H")) != string::npos)
	   {
		   str_act.allow_H();
		   changed = true;
		   log << " Allowing changing H " << endl;
	   }
	   if( simtype.find(string("ST")) != string::npos)
	   {
		   str_act.allow_ST();
		   changed = true;
		   log << " Allowing changing stacking " << endl;
	   }
	   if( simtype.find(string("CROSS")) != string::npos)
	   {
	  		str_act.allow_CROSS();
	  		changed = true;
	  		log << " Allowing changing crossstacking " << endl;
	   }
	   if( simtype.find(string("TEMP")) != string::npos)
	   {
	      str_act.allow_ST_dep();
	      changed = true;
	      log << " Allowing changing t-dep part of stacking " << endl;
	   }
	   if( simtype.find(string("AVG")) != string::npos)
	   {
	  	  average_only = true;
	  	  changed = true;
		  log << " Allowing changing only average! " << endl;
	   }


	   if(!changed)
		   throw string( string("Error processing simulation type, unrecognized option: ") + simtype);
	   str_new = str_best = str_act;
}
//----------------------------------------------------------------------------------
 Anneal_Simulation::Anneal_Simulation(string& parameter_file, const char *simtype, double _sim_beta,  double pmax, const char * logfile_prefix,  ostream& log )
: log(log),  simtype(simtype){

	   id = 0;
	   proc_no = 1;
	   string infoname =  string(logfile_prefix)  + string(simtype) + string(".log");
	   simulation_info_out.open(infoname.c_str());

	   string parname =  string(logfile_prefix)  + string(simtype) +  string(".par");
	   best_params_out.open(parname.c_str());
	   best_params_changed = false;

	   sim_beta = _sim_beta;
	   olddifference = 100;
	   newdifference = 100;
	   bestdifference = 100;
	   perturb_max = pmax;
	   int total_sequences = 0;

	//   correct_melts = 0;
	   correct_melts.clear();
	   attempt_melts = 0;
	   current_melts = 0;

	   log << "#Loading parameters from " << parameter_file << endl;
	   str_act.load_from_file(parameter_file.c_str());
	   cout << str_act << endl;
	   log << "#loaded: " << str_act << endl;

	   set_sim_type(string(simtype));
	   str_new = str_act;


}
//------------------------------------------------------------------------------------
void Anneal_Simulation::load_sequences_only(vector<string> &sequence_files, vector<string> &data_files)
{
	   int total_sequences = 0;
	   sequences_in_files.clear();

	   if(sequence_files.size() != data_files.size())
	   {
	     throw string("Incompatible dimensions ");
	   }

        if(id == 0)
	      log << "There are in total " << sequence_files.size() << endl;
	   for(int i = 0; i < sequence_files.size(); i++)
	   {
		 if(id == 0)
	         log << "Loading sequences  and energy data from " << sequence_files[i] << " and " << data_files[i] << endl;
	     int loaded = load_sequences_without_states(sequence_files[i].c_str(),data_files[i].c_str());
	     total_sequences += loaded;
	     for(int k = 0; k < loaded; k++)
	         sequences_in_files.push_back(i);

	   }

	   if(id == 0)
	     log << "Loaded number of sequences: " << total_sequences << endl;

	   //correct_melts = new double[total_sequences];
	   attempt_melts = new double[total_sequences];
	   current_melts = new double[total_sequences];

	   for(int i = 0; i < total_sequences; i++)
	   {
	         correct_melts[i] = states[i].get_correct_beta_melt(); //GetBetaForSequence(sequences[i].c_str(),states[i].get_box_size(),is_complementary_sequence(sequences[i]));
	         //str_act.set_relevant_parameters(sequences[i].c_str());
	   }


	   set_sim_type(string(simtype));

	   str_new = str_best = str_act;
}
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
int Anneal_Simulation::load_sequences_without_states(const char *seqfilename, const char *energyfilename)
{


   PAR_struct original_parameters;
   //ifstream seqfin(seqfilename);
   ifstream fin(energyfilename);

   if (!fin.good())
   {
     throw string ("Error loading data from ")+ string(energyfilename);
   }

   double beta_start, beta_step, beta0;
   int beta_array_length;
   double box;
   string simtype;
   fin >> simtype >> box >> beta0 >> beta_start >> beta_step >> beta_array_length;
   fin >> original_parameters;

   if(!fin.good())
   {
     throw(string( "Unable to load data from ") + string( energyfilename));
   }
   else
   {
     if(id == 0)
     {
	  log << "Loaded following simulation parameters from " << energyfilename << endl;
      log << " Box side: " << box << ", simulation run at beta=" << beta0 << " and will be fitted to " << beta_start << " up to " << beta_array_length << " * "  << beta_step << endl;
      log << " The loaded simulation was generated with parameters: " << original_parameters << endl;
     }
   }



   int start_from = states.size();
   int loaded_seq = load_sequences(seqfilename,beta0,beta_start,beta_step,beta_array_length,original_parameters,box,simtype);
   if(id == 0)
   {
	   log << "Loaded " << loaded_seq << " sequences from file " << seqfilename << endl;
	   simulation_info_out << "#Loaded data from  " << energyfilename << " and " << seqfilename << endl;

   }
   //now we load data into states structure


 return loaded_seq;

}
//-------------------------------------------------------------------------------------
void Anneal_Simulation::calculate_melts_one_sequence_at_time(vector<string> &data_files)
{
 //PAR_struct original_parameters;


 for(int i = 0; i < states.size(); i++)
 {
	 //states[i].clear();
	 int states_files_index = sequences_in_files[i];
     ifstream fin(data_files[states_files_index].c_str());
	 char buffer[4096];
     fin.getline(buffer,4096);
     fin.getline(buffer,4096); //parameters
     fin.getline(buffer,4096); //actual data
	 while(fin.good())
	 {
	     string str(buffer);
	     if(str.length() > 1 && str[0] != '\n')
	     {
	        states[i].add_state(str);
	        //cout << "Adding " << str << endl;
	     }
	      fin.getline(buffer,4096);
	 }
     fin.close();
     //cout << "loaded " << states[i].states.size() << " for " << i << " th sequences " << states[i].sequence << endl;


    // cout << str_act << endl;
	 attempt_melts[i] = current_melts[i] = states[i].calculate_beta_melt_for_parameters(str_act);

     states[i].clear();
 }
}
//-----------------------------------------------------------------------------------
//This loads only fraction of states calculated from myid and total_id (useful for MPI program)
//-------------------------------------------------------------------------------------
void Anneal_Simulation::load_fraction_of_states(vector<string> &data_files,int myid, int total_id)
{


 int states_to_process = states.size() / total_id;
 int stop = (myid+1)*states_to_process;
 if(myid == total_id - 1)
	  stop = states.size();
 for(int i = myid*states_to_process; i < stop; i++)
 {

	 int states_files_index = sequences_in_files[i];
     ifstream fin(data_files[states_files_index].c_str());
	 char buffer[4096];
     fin.getline(buffer,4096); //init line
     fin.getline(buffer,4096); //parameters
     fin.getline(buffer,4096); //actual data
	 while(fin.good())
	 {
	     string str(buffer);
	    // cout << " loaded " << str << endl;
	     if(str.length() > 1 && str[0] != '\n')
	     {
	        states[i].add_state(str);

	     }
	     fin.getline(buffer,4096);
	 }
     fin.close();
     //cout << "This is " << id << " and I loaded " << states[i].states.size() << " states into " << i << " th sequence " << states[i].sequence << endl;
	 attempt_melts[i] = current_melts[i] = states[i].get_correct_beta_melt();
	 //states[i].clear();
 }
}
//-------------------------------------------------------------------------------------
void Anneal_Simulation::calculate_attempt_melts(PAR_struct &pars, int myid, int total_id)
{
	int states_to_process = states.size() / total_id;
	int stop = (myid+1)*states_to_process;
	if(myid == total_id - 1)
		  stop = states.size();
	for(int i = myid*states_to_process; i < stop; i++)
	{
		  attempt_melts[i] = states[i].calculate_beta_melt_for_parameters(pars);
	}

}
//-------------------------------------------------------------------------------------
int Anneal_Simulation::calculate_attempt_melts(PAR_struct &pars)
 {
   for(int i = 0; i < states.size(); i++)
   {
     attempt_melts[i] = states[i].calculate_beta_melt_for_parameters(pars);
     if(attempt_melts[i] == 0)
    	 return 0;
   }
   return 1;
 }
//-------------------------------------------------------------------------------------
void Anneal_Simulation::calculate_attempt_yields(PAR_struct &pars, int myid, int total_id)
{
	int states_to_process = states.size() / total_id;
	int stop = (myid+1)*states_to_process;
	if(myid == total_id - 1)
		  stop = states.size();
	for(int i = myid*states_to_process; i < stop; i++)
	{
		  attempt_melts[i] = states[i].get_distance_from_melting_curve(pars);
	}

}
//-------------------------------------------------------------------------------------
void Anneal_Simulation::calculate_attempt_widths(PAR_struct &pars,int myid, int total_id)
{
	int states_to_process = states.size() / total_id;
	int stop = (myid+1)*states_to_process;
	if(myid == total_id - 1)
		stop = states.size();
	for(int i = myid*states_to_process; i < stop; i++)
	{
		attempt_melts[i] = states[i].get_width_of_melting_curve(pars);
	}

}
//------------------------------------------------------------------------------------
void Anneal_Simulation::show_yields(void)
{
 for(int i = 0; i < states.size(); i++)
 {
	 states[i].show_yields_for_parameters(str_act);
 }
}

//-------------------------------------------------------------------------------------
void Anneal_Simulation::calculate_attempt_yields(PAR_struct &pars)
 {
   for(int i = 0; i < states.size(); i++)
   {
     attempt_melts[i] = states[i].get_distance_from_melting_curve(pars);
   }

 }
//-------------------------------------------------------------------------------------
void Anneal_Simulation::calculate_attempt_widths(PAR_struct &pars)
{
	//cout << "Calculating attempt widrhs!" << endl;
   for(int i = 0; i < states.size(); i++)
   {
     attempt_melts[i] = states[i].get_width_of_melting_curve(pars);
   }

}
//-------------------------------------------------------------------------------------
void Anneal_Simulation::calculate_correct_widths(void)
{
   for(int i = 0; i < states.size(); i++)
   {
     correct_melts[i] = states[i].get_correct_width();
   }
}

//-------------------------------------------------------------------------------
void Anneal_Simulation::show_widths(ostream& out)
{
 //  calculate_attempt_widths(str_act);
 //  calculate_correct_widths();
  // if(total_sequences != states.size())
//	   throw string("erorr, the number of sequences doses not match the number of states");
   for(int i = 0; i < states.size(); i++)
   {
	 if(attempt_melts[i] < 0)
		 out << "#!!!Error, calculated attempt width is < 0 " << endl;
	 if(correct_melts[i] < 0)
	 		 out << "#!!!Error, calculated correct width is < 0 " << endl;

     out << sequences[i] << " " << attempt_melts[i] - correct_melts[i]  << " " << attempt_melts[i] << " " << correct_melts[i] << " " << states[i].get_correct_width() << endl;
   }

}
//-----------------------------------------------------------------------------------
double Anneal_Simulation::sum_differences(void)
 {
   double difference = 0;

   for(int i = 0; i < states.size(); i++)
   {

     difference += fabs(correct_melts[i] - attempt_melts[i] );// SQR(correct_melts[i] - attempt_melts[i] );
   }
   return difference / states.size();
 }
//-----------------------------------------------------------------------------------
double Anneal_Simulation::pow4_differences(void)
 {
   double difference = 0;

   for(int i = 0; i < states.size(); i++)
   {

     difference += (correct_melts[i] - attempt_melts[i] )*(correct_melts[i] - attempt_melts[i] ) * (correct_melts[i] - attempt_melts[i] )*(correct_melts[i] - attempt_melts[i] );// SQR(correct_melts[i] - attempt_melts[i] );
   }
   return difference / states.size();
 }

//-----------------------------------------------------------------------------------
int Anneal_Simulation::load_data_from_files(const char *seqfilename, const char *energyfilename)
 {
   
  
   PAR_struct original_parameters;
   //ifstream seqfin(seqfilename);
   ifstream fin(energyfilename);
   
   if (!fin.good())
   {
     throw string ("Error loading data from ")+ string(energyfilename);
   }

   double beta_start, beta_step, beta0;
   int beta_array_length;
   double box;
   string simtype;
   fin >> simtype >> box >> beta0 >> beta_start >> beta_step >> beta_array_length; 

#ifdef DBG
   cout << "Loaded simtype: " <<  simtype << " "<< box << " " << beta0 << " " << beta_start << " " << beta_step << " " <<  " " << beta_array_length << endl;
#endif

   fin >> original_parameters;


   if(!fin.good())
   {
     throw(string( "Unable to load data from ") + string( energyfilename));
   }
   else
   {
     log << "#Loaded following simulation parameters from " << energyfilename << endl;
     log << "# Box side: " << box << ", simulation run at beta=" << beta0 << " and will be fitted to " << beta_start << " up to " << beta_array_length << " * "  << beta_step << endl;
     log << "# The loaded simulation was generated with parameters: " << original_parameters << endl;
   }

   

   int start_from = states.size();
   int loaded_seq = load_sequences(seqfilename,beta0,beta_start,beta_step,beta_array_length,original_parameters,box,simtype);
    log << "#Loaded " << loaded_seq << " sequences from file " << seqfilename << endl;
   //now we load data into states structure
   log << "#Loading all simulation states  from " << energyfilename << endl;
   char buffer[4096];
   fin.getline(buffer,4096);
   while(fin.good())
   {
    string str(buffer);
    if(str.length() > 1 && str[0] != '\n') 
    {
     for(int i = start_from; i < start_from + loaded_seq ; i++)
     {
     
       //cout << "Adding " << str << endl;

       states[i].add_state(str);
       
     }
    }
     fin.getline(buffer,4096);   
   }

  simulation_info_out << "#Loaded data from  " << energyfilename << " and " << seqfilename << " and started with parameters " << original_parameters << endl;

  log << " simulation states loaded " << endl;
  
 return loaded_seq;

 }
//--------------------------------------------------------------------------------
double Anneal_Simulation::accept_move(void)
 {
  for(int i = 0; i < states.size(); i++)
  {
    current_melts[i] = attempt_melts[i];
  }
  str_act = str_new;
  olddifference = newdifference;  

  if(newdifference < bestdifference)
  {
    bestdifference = newdifference;
    str_best = str_new;
  }

  return newdifference;
 }
 //-----------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------
 int Anneal_Simulation::simulation_STUN_MC_step(double sim_beta)
 {
   int accepted = 0;
   int notall_in_limits = 0;
   do
   {
     str_new.perturb_parameter(perturb_max);
     notall_in_limits = 0;
     for(int i = 0; i < 2; i++)
     {
    	 if(str_new.vals[i] < 0.7 || str_new.vals[i] > 2.2)
    	       {
    	    	   notall_in_limits = 1;
    	    	   //break;
    	       }
     }
     for(int i =2; i < PAR_SIZE; i++)
     {
       if(str_new.vals[i] < 0.8 || str_new.vals[i] > 3.5)
       {
    	   notall_in_limits = 1;
    	   //break;
       }
     }
     if(notall_in_limits)
    	 str_new = str_act;

   } while(notall_in_limits);
  // cout << "Parspert " << str_new << endl;
   double _gama = 10.;
   calculate_attempt_melts(str_new);
   newdifference = sum_differences();

   double e_stun_new = 1.0 - exp(-_gama * (newdifference - bestdifference) );
   double e_stun_old = 1.0 - exp(-_gama * (olddifference - bestdifference) );

   if ( newdifference < olddifference || rand_01()  < exp(-sim_beta * (e_stun_new - e_stun_old))   )
   {
     accept_move();
     accepted = 1;

   }
   else
   {
     str_new = str_act;
     accepted = 0;
   }

  // cout << "Parspert " << accepted << "  " << newdifference <<  " " << olddifference << endl;
   return accepted;
 }
//-------------------------------------------------------------------------------------
 int Anneal_Simulation::simulation_MC_step(double sim_beta)
 {
   int accepted = 0;
   str_new.perturb_parameter(perturb_max,average_only);

   if(!str_new.check_reasonable_limits())
   {
	   accepted = 0;
	   str_new = str_act;
	   return 0;
   }
  // cout << "Parspert " << str_new << endl;

   calculate_attempt_melts(str_new);
   newdifference = sum_differences();
   if ( newdifference < olddifference || rand_01()  < exp(-sim_beta * (newdifference - olddifference))   )
   {
     accept_move();
     accepted = 1;

   }  
   else
   {
     str_new = str_act;
     accepted = 0;
   }

  // cout << "Parspert " << accepted << "  " << newdifference <<  " " << olddifference << endl;
   return accepted;
 } 

 //-------------------------------------------------------------------------------------
 //-------------------------------------------------------------------------------------
  int Anneal_Simulation::simulation_MC_step_average(double sim_beta)
  {
    int accepted = 0;


    str_new.perturb_parameter_average(perturb_max);
    if(!str_new.check_reasonable_limits())
    {
    	accepted = 0;
    	return accepted;
    }
   // cout << "Parspert " << str_new << endl;

    if( ! calculate_attempt_melts(str_new) )
    {
    	accepted = 0;
    	return 0;
    }

    newdifference = sum_differences();
    if ( newdifference < olddifference || rand_01()  < exp(-sim_beta * (newdifference - olddifference))   )
    {
      accept_move();
      accepted = 1;

    }
    else
    {
      str_new = str_act;
      accepted = 0;
    }

   // cout << "Parspert " << accepted << "  " << newdifference <<  " " << olddifference << endl;
    return accepted;
  }

  //-------------------------------------------------------------------------------------

 //-------------------------------------------------------------------------------------
  int Anneal_Simulation::simulation_MC_step_with_yields(double sim_beta)
  {
    int accepted = 0;
    str_new.perturb_parameter(perturb_max);

   // cout << "Parspert " << str_new << endl;

    calculate_attempt_yields(str_new);
    newdifference = sum_differences();
    if ( newdifference < olddifference || rand_01()  < exp(-sim_beta * (newdifference - olddifference))   )
    {
      accept_move();
      accepted = 1;

    }
    else
    {
      str_new = str_act;
      accepted = 0;
    }

   // cout << "Parspert " << accepted << "  " << newdifference <<  " " << olddifference << endl;
    return accepted;
  }
//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
  int Anneal_Simulation::simulation_MC_step_with_widths(double sim_beta)
  {
	  int accepted = 0;
	  str_new.perturb_parameter(perturb_max,average_only);

	  // cout << "Parspert " << str_new << endl;

	  calculate_attempt_widths(str_new);
	  newdifference = sum_differences();
	  if ( newdifference < olddifference || rand_01()  < exp(-sim_beta * (newdifference - olddifference))   )
	  {
		  accept_move();
		  accepted = 1;

	  }
	  else
	  {
		  str_new = str_act;
		  accepted = 0;
	  }

	  // cout << "Parspert " << accepted << "  " << newdifference <<  " " << olddifference << endl;
	  return accepted;
  }
 //-------------------------------------------------------------------------------------

void Anneal_Simulation::save_best_results(int iter)
 { 
   best_params_out << "#Saving best result achieved by " << iter << " iterations "<< " with Hamiltonian " << bestdifference <<  endl;
   str_best.show_parameters(best_params_out);
 }


//---------------------------------------------------------------------------------
void Anneal_Simulation::simulated_annealing(int steps, double beta_start, double high_acceptance , double low_acceptance , int adapt_after, int output_log_after )
 {
   double pert;
   int accepted_moves = 0;
   double beta = beta_start;
   int accept_one = 0;
   double maxdif = 0;
   double avgdif = 0;

   calculate_attempt_melts(str_new);
   olddifference = bestdifference = newdifference = sum_differences();

   maxdif = get_max_difference_kelvin();
   avgdif =  get_avg_difference_kelvin();

   log << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
   simulation_info_out << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
   for(int iter = 1; iter <= steps; iter++)
   {
     accept_one = simulation_MC_step(beta);
     accepted_moves += accept_one;
     if(accept_one)
     {
    	 maxdif = get_max_difference_kelvin();
    	 avgdif =  get_avg_difference_kelvin();
     }
     log << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif<< " " << maxdif << " " << str_act << endl;
  
     simulation_info_out << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif << " " << maxdif << " " << str_act << endl;
   

     if(iter % adapt_after == 0)
     {
       if (double(accepted_moves) / adapt_after  < low_acceptance)
       {
          pert = rand_01() * 0.1;  
          beta /= 1. + pert;
       }

       if (double(accepted_moves) / adapt_after  > high_acceptance)
       {
          pert = rand_01() * 0.1;  
          beta *= 1. + pert;
       }        
       accepted_moves = 0;

     }     
    
     if(iter % output_log_after == 0)
     {
       log << "Saving best results so far " << endl;
       save_best_results(iter);
     }
   }

   
  
 }
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
void Anneal_Simulation::simulated_annealing_average(int steps, double beta_start, double high_acceptance , double low_acceptance , int adapt_after, int output_log_after )
 {
   double pert;
   int accepted_moves = 0;
   double beta = beta_start;
   int accept_one = 0;
   double maxdif = 0;
   double avgdif = 0;

   calculate_attempt_melts(str_new);
   olddifference = bestdifference = newdifference = sum_differences();

   maxdif = get_max_difference_kelvin();
   avgdif =  get_avg_difference_kelvin();

   log << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
   simulation_info_out << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
   for(int iter = 1; iter <= steps; iter++)
   {
     accept_one = simulation_MC_step_average(beta);
     accepted_moves += accept_one;
     if(accept_one)
     {
    	 maxdif = get_max_difference_kelvin();
    	 avgdif =  get_avg_difference_kelvin();
     }
     log << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif<< " " << maxdif << " " << str_act << endl;

     simulation_info_out << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif << " " << maxdif << " " << str_act << endl;


     if(iter % adapt_after == 0)
     {
       if (double(accepted_moves) / adapt_after  < low_acceptance)
       {
          pert = rand_01() * 0.1;
          beta /= 1. + pert;
       }

       if (double(accepted_moves) / adapt_after  > high_acceptance)
       {
          pert = rand_01() * 0.1;
          beta *= 1. + pert;
       }
       accepted_moves = 0;

     }

     if(iter % output_log_after == 0)
     {
       log << "Saving best results so far " << endl;
       save_best_results(iter);
     }
   }



 }
//---------------------------------------------------------------------------------
void Anneal_Simulation::simulated_annealing_with_widths(int steps, double beta_start, double high_acceptance , double low_acceptance , int adapt_after, int output_log_after )
 {
   double pert;
   int accepted_moves = 0;
   double beta = beta_start;
   int accept_one = 0;
   double maxdif = 0;
   double avgdif = 0;

   //calculate_attempt_melts(str_new);
   calculate_attempt_widths(str_new);
   for(int i = 0; i < states.size(); i++)
   {
	  // correct_melts[i] = states[i].get_correct_width();
	   current_melts[i] = attempt_melts[i];
        //   cout << "For seq " << i << " the distance is " << attempt_melts[i] << endl;
	   //states[i].fill_in_predicted_yields();
   }
   olddifference = bestdifference = newdifference = sum_differences();

   maxdif = get_max_difference();
   avgdif =  get_avg_difference();

   log << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
   simulation_info_out << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
   for(int iter = 1; iter <= steps; iter++)
   {
     accept_one = simulation_MC_step_with_widths(beta);
     accepted_moves += accept_one;
     if(accept_one)
     {
    	 maxdif = get_max_difference();
    	 avgdif =  get_avg_difference();
     }
     log << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif<< " " << maxdif << " " << str_act << endl;

     simulation_info_out << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif << " " << maxdif << " " << str_act << endl;


     if(iter % adapt_after == 0)
     {
       if (double(accepted_moves) / adapt_after  < low_acceptance)
       {
          pert = rand_01() * 0.1;
          beta /= 1. + pert;
       }

       if (double(accepted_moves) / adapt_after  > high_acceptance)
       {
          pert = rand_01() * 0.1;
          beta *= 1. + pert;
       }
       accepted_moves = 0;

     }

     if(iter % output_log_after == 0)
     {
       log << "Saving best results so far " << endl;
       save_best_results(iter);
     }
   }



 }
//---------------------------------------------------------------------------------
void Anneal_Simulation::simulated_annealing_with_yields(int steps, double beta_start, double high_acceptance , double low_acceptance , int adapt_after, int output_log_after )
 {
   double pert;
   int accepted_moves = 0;
   double beta = beta_start;
   int accept_one = 0;
   double maxdif = 0;
   double avgdif = 0;

   //calculate_attempt_melts(str_new);
   calculate_attempt_yields(str_new);
   for(int i = 0; i < states.size(); i++)
   {
	   correct_melts[i] = 0;
	   current_melts[i] = attempt_melts[i];
        //   cout << "For seq " << i << " the distance is " << attempt_melts[i] << endl;
	   //states[i].fill_in_predicted_yields();
   }
   olddifference = bestdifference = newdifference = sum_differences();

   maxdif = get_max_difference();
   avgdif =  get_avg_difference();

   log << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
   simulation_info_out << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
   for(int iter = 1; iter <= steps; iter++)
   {
     accept_one = simulation_MC_step_with_yields(beta);
     accepted_moves += accept_one;
     if(accept_one)
     {
    	 maxdif = get_max_difference();
    	 avgdif =  get_avg_difference();
     }
     log << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif<< " " << maxdif << " " << str_act << endl;

     simulation_info_out << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif << " " << maxdif << " " << str_act << endl;


     if(iter % adapt_after == 0)
     {
       if (double(accepted_moves) / adapt_after  < low_acceptance)
       {
          pert = rand_01() * 0.1;
          beta /= 1. + pert;
       }

       if (double(accepted_moves) / adapt_after  > high_acceptance)
       {
          pert = rand_01() * 0.1;
          beta *= 1. + pert;
       }
       accepted_moves = 0;

     }

     if(iter % output_log_after == 0)
     {
       log << "Saving best results so far " << endl;
       save_best_results(iter);
     }
   }



 }
//----------------------------------------------------------------------------------
void Anneal_Simulation::STUN_annealing(int steps, double beta_start, double high_acceptance, double low_acceptance , int adapt_after, int output_log_after)
 {
   double pert;
   int accepted_moves = 0;
   double beta = beta_start;
   int accept_one = 0;
   double maxdif = 0;
   double avgdif = 0;
   int adaptations = 0;
 //  double change = double(high_temp - low_temp) / adapt_step;

   calculate_attempt_melts(str_new);
   olddifference = bestdifference = newdifference = sum_differences();

   maxdif = get_max_difference_kelvin();
   avgdif =  get_avg_difference_kelvin();

   log << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
   simulation_info_out << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
   for(int iter = 1; iter <= steps; iter++)
   {
     accept_one = simulation_STUN_MC_step(beta);
     accepted_moves += accept_one;
     if(accept_one)
     {
    	 maxdif = get_max_difference_kelvin();
    	 avgdif =  get_avg_difference_kelvin();
     }
     log << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif<< " " << maxdif << " " << str_act << endl;

     simulation_info_out << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif << " " << maxdif << " " << str_act << endl;


     if(iter % adapt_after == 0)
     {
    	 if (double(accepted_moves) / adapt_after  < low_acceptance)
    	 {
    		 pert = rand_01() * 0.1;
    		 beta /= 1. + pert;
    	 }

    	 if (double(accepted_moves) / adapt_after  > high_acceptance)
    	 {
    		 pert = rand_01() * 0.1;
    		 beta *= 1. + pert;
    	 }
    	 accepted_moves = 0;
     }

     if(iter % output_log_after == 0)
     {
       log << "Saving best results so far " << endl;
       save_best_results(iter);
     }
   }



 }
//---------------------------------------------------------------------------------
void Anneal_Simulation::retherm_simulated_annealing(int steps, double beta_start, double high_temp, double low_temp , int adapt_after,int adapt_step, int output_log_after)
 {
   double pert;
   int accepted_moves = 0;
   double beta = beta_start/high_temp;
   int accept_one = 0;
   double maxdif = 0;
   double avgdif = 0;
   int adaptations = 0;
   double change = double(high_temp - low_temp) / adapt_step;

   calculate_attempt_melts(str_new);
   olddifference = bestdifference = newdifference = sum_differences();

   maxdif = get_max_difference_kelvin();
   avgdif =  get_avg_difference_kelvin();

   log << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
   simulation_info_out << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
   for(int iter = 1; iter <= steps; iter++)
   {
     accept_one = simulation_MC_step(beta);
     accepted_moves += accept_one;
     if(accept_one)
     {
    	 maxdif = get_max_difference_kelvin();
    	 avgdif =  get_avg_difference_kelvin();
     }
     log << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif<< " " << maxdif << " " << str_act << endl;

     simulation_info_out << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif << " " << maxdif << " " << str_act << endl;


     if(iter % adapt_after == 0)
     {
       adaptations++;
       if(adaptations >= adapt_step)
       {
    	   adaptations = 0;
    	   beta = beta_start / high_temp;
       }
       else
       {
    	   beta = beta_start / (high_temp - adaptations*change);
       }

       accepted_moves = 0;

     }

     if(iter % output_log_after == 0)
     {
       log << "Saving best results so far " << endl;
       save_best_results(iter);
     }
   }



 }
//----------------------------------------------------------------------------------
void Anneal_Simulation::show_all_melts(ostream& out, double cutoff)
 {
   out << "#Sequence (difference in kelvin)  melt correct_melt " << endl;
   double diff;
   for(int i =0; i < states.size(); i++)
   {
     diff = beta_to_kelvin(current_melts[i]) - beta_to_kelvin(correct_melts[i]);
     if(cutoff <= 0 || fabs(diff) > cutoff) 
        out << sequences[i] << " " <<   diff << " " << beta_to_kelvin(current_melts[i]) << " " << beta_to_kelvin(correct_melts[i]) << endl;
     
   }

   out << "#Maxdifference is: " << get_max_difference_kelvin() << " and average difference is: " << get_avg_difference_kelvin() << endl;

 }

//--------------------------------------------------------------------------------
ostream &Anneal_Simulation::show(ostream& out)
{
   for(int i = 0; i < states.size(); i++)
   {
     out << states[i] << endl;
   }
  
   return out;

}

//---------------------------------------------------------------------------------

 double Anneal_Simulation::get_max_difference(void)
 {
   double maxdifference = 0;
   double difr;

   for(int i = 0; i < states.size(); i++)
   {
     difr = fabs(correct_melts[i] - attempt_melts[i] );
     if ( difr > maxdifference ) 
             maxdifference = difr  ;
   }
   return maxdifference;
 }
//---------------------------------------------------------------------------------
 double Anneal_Simulation::get_max_difference_kelvin(void)
 {
   double maxdifference = 0;
   double difr;

   for(int i = 0; i < states.size(); i++)
   {
     difr = fabs(beta_to_kelvin(correct_melts[i]) - beta_to_kelvin(attempt_melts[i]) );
     if ( difr > maxdifference ) 
             maxdifference = difr  ;
   }
   return maxdifference;
 }

 //---------------------------------------------------------------------------------
 double Anneal_Simulation::get_avg_difference_kelvin(void)
 {
   double avgdifference = 0;
   double difr;

   for(int i = 0; i < states.size(); i++)
   {
     avgdifference += fabs(beta_to_kelvin(correct_melts[i]) - beta_to_kelvin(attempt_melts[i]) );
   }
  
   return avgdifference /= states.size();

 }
 //---------------------------------------------------------------------------------
  double Anneal_Simulation::get_avg_difference(void)
  {
    double avgdifference = 0;
    double difr;

    for(int i = 0; i < states.size(); i++)
    {
      avgdifference += fabs((correct_melts[i]) - (attempt_melts[i]) );
    }

    return avgdifference /= states.size();

  }
//---------------------------------------------------------------------------------
void Anneal_Simulation::simulated_annealing_with_maxdifference(int steps, double beta_start, double high_acceptance , double low_acceptance , int adapt_after, int output_log_after )
 {
   double pert;
   int accepted_moves = 0;
   double beta = beta_start;
   olddifference = bestdifference = newdifference = get_max_difference();

   log << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
   simulation_info_out << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
   for(int iter = 1; iter <= steps; iter++)
   {
     accepted_moves += simulation_MC_step_with_maxdifference(beta);
    
    log << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << get_avg_difference_kelvin() << " " << get_max_difference_kelvin() << " " << str_act << endl;
  
     simulation_info_out << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << get_avg_difference_kelvin() << " " << get_max_difference_kelvin() << " " << str_act << endl;
   

     if(iter % adapt_after == 0)
     {
       if (double(accepted_moves) / adapt_after  < low_acceptance)
       {
          pert = rand_01() * 0.1;  
          beta /= 1. + pert;
       }

       if (double(accepted_moves) / adapt_after  > high_acceptance)
       {
          pert = rand_01() * 0.1;  
          beta *= 1. + pert;
       }        
       accepted_moves = 0;

     }     
    
     if(iter % output_log_after == 0)
     {
       log << "Saving best results so far " << endl;
       save_best_results(iter);
     }
   }

   
  
 }
//----------------------------------------------------------------------------------

 int Anneal_Simulation::simulation_MC_step_with_maxdifference(double sim_beta)
 {
   int accepted = 0;
   str_new.perturb_parameter(perturb_max);
   calculate_attempt_melts(str_new);
   newdifference = get_max_difference();
   if ( newdifference < olddifference || rand_01()  < exp(-sim_beta * (newdifference - olddifference))   )
   {
     accept_move();
     accepted = 1;
   }  
   else
   {
     str_new = str_act;
     accepted = 0;
   }

   return accepted;
 } 
 //-------------------------------------------------------------------------------------

 /*
 //----------------------------------------------------------------------------------
  int Anneal_Simulation::simulation_MC_step_with_avg_ratios(double sim_beta)
  {
    int accepted = 0;
    str_new.perturb_avg_ratio(perturb_max);
    calculate_attempt_melts(str_new);
    newdifference = sum_differences();
    if ( newdifference < olddifference || rand_01()  < exp(-sim_beta * (newdifference - olddifference))   )
    {
      accept_move();
      accepted = 1;
    }
    else
    {
      str_new = str_act;
      accepted = 0;
    }

    return accepted;
  }
  */

//-------------------------------------------------------------------------------------
#ifdef HAVE_MPI
 Anneal_Simulation::Anneal_Simulation(int rank, int total_proc,vector<string> &sequence_files, vector<string> &data_files,string& parameter_file, const char *simtype, double _sim_beta,  double pmax, const char * logfile_prefix,  ostream& log )
: log(log),  simtype(simtype) {
       id = rank;
       proc_no = total_proc;
       cout << "This is process " << id << endl;

       if(id == 0)
       {
 	    string infoname =  string(logfile_prefix)  + string(simtype) + string(".log");
 	    simulation_info_out.open(infoname.c_str());

 	    string parname =  string(logfile_prefix)  + string(simtype) +  string(".par");
 	    best_params_out.open(parname.c_str());
 	    best_params_changed = false;
       }
 	   sim_beta = _sim_beta;
 	   olddifference = 100;
 	   newdifference = 100;
 	   bestdifference = 100;
 	   perturb_max = pmax;

 	   int total_sequences = 0;

 	   //correct_melts = 0;
 	   attempt_melts = 0;
 	   current_melts = 0;

 	  if(id == 0)
      {
 	   log << "Loading parameters from " << parameter_file << endl;
 	   str_act.load_from_file(parameter_file.c_str());
 	   log << "loaded " << str_act << endl;
 	   simulation_info_out  << "# Started with parameters " << str_act << endl;
#ifdef NO_T_DEP
       simulation_info_out << "#No T dependence in stacking!!!" << endl;
#endif
      }



 	 load_sequences_only(sequence_files,data_files);
 	 load_fraction_of_states(data_files,id,proc_no);


 }
 //---------------------------------------------------------------------------------------
int Anneal_Simulation::MPI_simulation_MC_step(double sim_beta)
{
  int accepted = 0;
  if(id == 0)
  {
	str_new.perturb_parameter(perturb_max);
	//cout << "Parspert " << str_new << endl;
	MPI_send_pars_new();
  }
  else
  {
    MPI_receive_pars_new();
  }

  calculate_attempt_melts(str_new,id,proc_no);
  if(id == 0)
  {
	  MPI_receive_melts();
  }
  else
  {
	  MPI_send_melts();
  }
  if(id == 0)
  {
	   newdifference = sum_differences();
	   if ( newdifference < olddifference || rand_01()  < exp(-sim_beta * (newdifference - olddifference))   )
	   {
	     accept_move();
	     accepted = 1;
	   }
	   else
	   {
	     str_new = str_act;
	     accepted = 0;
	   }
	  // cout << "Parspert " << accepted << "  " << newdifference << " " << olddifference <<  endl;
  }

  return accepted;

}
//-----------------------------------------------------------------------------------------
void Anneal_Simulation::MPI_simulated_annealing(int steps, double beta_start, double high_acceptance , double low_acceptance , int adapt_after, int output_log_after )
{
   double pert;
   int accepted_moves = 0;
   double beta = beta_start;
   int accept_one = 0;
   double maxdif = 0;
   double avgdif = 0;

   if(id == 0)
   {
   	  MPI_send_pars_new();
   }
   else
   {
       MPI_receive_pars_new();
   }

   calculate_attempt_melts(str_new,id,proc_no);

   if(id == 0)
   {
   	  MPI_receive_melts();
   }
   else
   {
   	  MPI_send_melts();
   }
   if(id == 0)
   {
	  for(int i =0; i < states.size(); i++)
		   current_melts[i] = attempt_melts[i];
      olddifference = bestdifference = newdifference = sum_differences();

      maxdif =  get_max_difference_kelvin();
      avgdif = get_avg_difference_kelvin();
      log << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
      simulation_info_out << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
   }

   for(int iter = 1; iter <= steps; iter++)
   {

	accept_one = MPI_simulation_MC_step(beta);
    accepted_moves += accept_one ;

    if(id == 0)
    {
      if(accept_one)
      {
    	maxdif =  get_max_difference_kelvin();
    	avgdif = get_avg_difference_kelvin();
      }

      log << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif << " " << maxdif << " " << str_act << endl;
      simulation_info_out << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif << " " << maxdif << " " << str_act << endl;


      if(iter % adapt_after == 0)
      {
       if (double(accepted_moves) / adapt_after  < low_acceptance)
       {
          pert = rand_01() * 0.1;
          beta /= 1. + pert;
       }

       if (double(accepted_moves) / adapt_after  > high_acceptance)
       {
          pert = rand_01() * 0.1;
          beta *= 1. + pert;
       }
       accepted_moves = 0;

      }

      if(iter % output_log_after == 0 )
      {
       log << "Saving best results so far " << endl;
       save_best_results(iter);
      }
    }
   }



}
//------------------------------------------------------------------------------------------------------------------
void Anneal_Simulation::MPI_simulated_annealing_with_avg_ratios(int steps, double beta_start, double high_acceptance , double low_acceptance , int adapt_after, int output_log_after )
{
   double pert;
   int accepted_moves = 0;
   double beta = beta_start;
   int accept_one = 0;
   double maxdif = 0;
   double avgdif = 0;

   if(id == 0)
   {
   	  MPI_send_pars_new();
   }
   else
   {
       MPI_receive_pars_new();
   }

   calculate_attempt_melts(str_new,id,proc_no);

   if(id == 0)
   {
   	  MPI_receive_melts();
   }
   else
   {
   	  MPI_send_melts();
   }
   if(id == 0)
   {
	  for(int i =0; i < states.size(); i++)
	  {
		   current_melts[i] = attempt_melts[i];
		   correct_melts[i] = states[i].get_avg_beta_melt();
	  }
      olddifference = bestdifference = newdifference = sum_differences();

      maxdif =  get_max_difference_kelvin();
      avgdif = get_avg_difference_kelvin();
      log << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
      simulation_info_out << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
   }

   for(int iter = 1; iter <= steps; iter++)
   {

    accept_one = MPI_simulation_MC_perturb_avg_ratio_step(beta);
    accepted_moves += accept_one ;

    if(id == 0)
    {
      if(accept_one)
      {
    	maxdif =  get_max_difference_kelvin();
    	avgdif = get_avg_difference_kelvin();
      }

      log << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif << " " << maxdif << " " << str_act << endl;
      simulation_info_out << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif << " " << maxdif << " " << str_act << endl;


      if(iter % adapt_after == 0)
      {
       if (double(accepted_moves) / adapt_after  < low_acceptance)
       {
          pert = rand_01() * 0.1;
          beta /= 1. + pert;
       }

       if (double(accepted_moves) / adapt_after  > high_acceptance)
       {
          pert = rand_01() * 0.1;
          beta *= 1. + pert;
       }
       accepted_moves = 0;

      }

      if(iter % output_log_after == 0 )
      {
       log << "Saving best results so far " << endl;
       save_best_results(iter);
      }
    }
   }



}
//-----------------------------------------------------------------------------------------
void Anneal_Simulation::MPI_STUN_annealing(int steps, double beta_start, double high_acceptance , double low_acceptance , int adapt_after, int output_log_after )
{
	double pert;
	   int accepted_moves = 0;
	   double beta = beta_start;
	   int accept_one = 0;
	   double maxdif = 0;
	   double avgdif = 0;

	   if(id == 0)
	   {
	   	  MPI_send_pars_new();
	   }
	   else
	   {
	       MPI_receive_pars_new();
	   }

	   calculate_attempt_melts(str_new,id,proc_no);

	   if(id == 0)
	   {
	   	  MPI_receive_melts();
	   }
	   else
	   {
	   	  MPI_send_melts();
	   }
	   if(id == 0)
	   {
		  for(int i =0; i < states.size(); i++)
			   current_melts[i] = attempt_melts[i];
	      olddifference = bestdifference = newdifference = sum_differences();

	      maxdif =  get_max_difference_kelvin();
	      avgdif = get_avg_difference_kelvin();
	      log << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
	      simulation_info_out << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
	   }

	   for(int iter = 1; iter <= steps; iter++)
	   {

		accept_one = MPI_simulation_STUN_MC_step(beta);
	    accepted_moves += accept_one ;

	    if(id == 0)
	    {
	      if(accept_one)
	      {
	    	maxdif =  get_max_difference_kelvin();
	    	avgdif = get_avg_difference_kelvin();
	      }

	      log << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif << " " << maxdif << " " << str_act << endl;
	      simulation_info_out << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif << " " << maxdif << " " << str_act << endl;


	      if(iter % adapt_after == 0)
	      {
	       if (double(accepted_moves) / adapt_after  < low_acceptance)
	       {
	          pert = rand_01() * 0.1;
	          beta /= 1. + pert;
	       }

	       if (double(accepted_moves) / adapt_after  > high_acceptance)
	       {
	          pert = rand_01() * 0.1;
	          beta *= 1. + pert;
	       }
	       accepted_moves = 0;

	      }

	      if(iter % output_log_after == 0 )
	      {
	       log << "Saving best results so far " << endl;
	       save_best_results(iter);
	      }
	    }
	   }



}
//-----------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Anneal_Simulation::MPI_retherm_annealing(int steps, double beta_start, double high_temp , double low_temp, int adapt_after, int adapt_step, int output_log_after)
{
   double pert;
   int accepted_moves = 0;
   double beta = beta_start/high_temp;
   //double betamax = beta_start/high_temp;
   int accept_one = 0;
   double maxdif = 0;
   double avgdif = 0;
   int adaptations = 0;
   double change = double(-beta_start/high_temp + beta_start/low_temp) / adapt_step;

   if(id == 0)
   {
   	  MPI_send_pars_new();
   }
   else
   {
       MPI_receive_pars_new();
   }

   calculate_attempt_melts(str_new,id,proc_no);

   if(id == 0)
   {
   	  MPI_receive_melts();
   }
   else
   {
   	  MPI_send_melts();
   }
   if(id == 0)
   {
	  for(int i =0; i < states.size(); i++)
		   current_melts[i] = attempt_melts[i];
      olddifference = bestdifference = newdifference = sum_differences();

      maxdif =  get_max_difference_kelvin();
      avgdif = get_avg_difference_kelvin();
      log << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
      simulation_info_out << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
   }

   for(int iter = 1; iter <= steps; iter++)
   {

	accept_one = MPI_simulation_MC_step(beta);
    accepted_moves += accept_one ;

    if(id == 0)
    {
      if(accept_one)
      {
    	maxdif =  get_max_difference_kelvin();
    	avgdif = get_avg_difference_kelvin();
      }

      log << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif << " " << maxdif << " " << str_act << endl;
      simulation_info_out << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif << " " << maxdif << " " << str_act << endl;


      if(iter % adapt_after == 0)
      {

    	  adaptations++;
    	  if(adaptations > adapt_step)
    	  {
    		  adaptations = 0;
    		  beta = beta_start / high_temp;
    	  }
    	  else
    	  {
    		  beta = beta + change;
    	  }

    	  accepted_moves = 0;

      }

      if(iter % output_log_after == 0 )
      {
       log << "Saving best results so far " << endl;
       save_best_results(iter);
      }
    }
   }



}
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Anneal_Simulation::MPI_simulated_annealing_with_maxdifference(int steps, double beta_start, double high_acceptance , double low_acceptance , int adapt_after, int output_log_after )
{
   double pert;
   int accepted_moves = 0;
   double beta = beta_start;
   int accept_one = 0;

   if(id == 0)
   {
   	  MPI_send_pars_new();
   }
   else
   {
       MPI_receive_pars_new();
   }

   calculate_attempt_melts(str_new,id,proc_no);

   double maxdif =  0;
   double avgdif = 0;

   if(id == 0)
   {
   	  MPI_receive_melts();
   }
   else
   {
   	  MPI_send_melts();
   }
   if(id == 0)
   {
	  for(int i =0; i < states.size(); i++)
		   current_melts[i] = attempt_melts[i];
      olddifference = bestdifference = newdifference = get_max_difference();

      log << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
      simulation_info_out << "#simbeta accratio simH AvgError MaxError  parameters " << endl;

      maxdif =  get_max_difference_kelvin();
      avgdif = get_avg_difference_kelvin();
   }

   for(int iter = 1; iter <= steps; iter++)
   {

    accept_one = MPI_simulation_MC_step_maxdifference(beta);
    accepted_moves += accept_one;
    if(id == 0)
    {
      if(accept_one)
      {
    	   maxdif =  get_max_difference_kelvin();
    	   avgdif = get_avg_difference_kelvin();
      }
      log << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif << " " << maxdif << " " << str_act << endl;
      simulation_info_out << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif << " " << maxdif << " " << str_act << endl;


      if(iter % adapt_after == 0)
      {
       if (double(accepted_moves) / adapt_after  < low_acceptance)
       {
          pert = rand_01() * 0.1;
          beta /= 1. + pert;
       }

       if (double(accepted_moves) / adapt_after  > high_acceptance)
       {
          pert = rand_01() * 0.1;
          beta *= 1. + pert;
       }
       accepted_moves = 0;

      }

      if(iter % output_log_after == 0 )
      {
       log << "Saving best results so far " << endl;
       save_best_results(iter);
      }
    }
   }



}
//-----------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
void Anneal_Simulation::MPI_simulated_annealing_with_yields(int steps, double beta_start, double high_acceptance , double low_acceptance , int adapt_after, int output_log_after )
 {
   double pert;
   int accepted_moves = 0;
   double beta = beta_start;
   int accept_one = 0;
   double maxdif = 0;
   double avgdif = 0;


   if(id == 0)
   {
   	  MPI_send_pars_new();
   }
   else
   {
       MPI_receive_pars_new();
   }

   //calculate_attempt_melts(str_new);
   calculate_attempt_yields(str_new,id,proc_no);

   if(id == 0)
   {
	   MPI_receive_melts();
   }
   else
   {
	   MPI_send_melts();
   }

   if(id == 0)
   {
	   for(int i = 0; i < states.size(); i++)
	   {
		   correct_melts[i] = 0;
		   current_melts[i] = attempt_melts[i];
		   //cout << "For seq " << i << " the distance is " << attempt_melts[i] << endl;
		   //states[i].fill_in_predicted_yields();
	   }
	   olddifference = bestdifference = newdifference = sum_differences();

	   maxdif = get_max_difference();
	   avgdif =  get_avg_difference();

	   log << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
	   simulation_info_out << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
   }

   for(int iter = 1; iter <= steps; iter++)
   {
	   accept_one = MPI_simulation_MC_step_with_yields(beta);
	   accepted_moves += accept_one;
	   if(id == 0)
	   {
		   if(accept_one)
		   {
			   maxdif = get_max_difference();
			   avgdif =  get_avg_difference();
		   }
		   log << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif<< " " << maxdif << " " << str_act << endl;

		   simulation_info_out << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif << " " << maxdif << " " << str_act << endl;


		   if(iter % adapt_after == 0)
		   {
			   if (double(accepted_moves) / adapt_after  < low_acceptance)
			   {
				   pert = rand_01() * 0.1;
				   beta /= 1. + pert;
			   }

			   if (double(accepted_moves) / adapt_after  > high_acceptance)
			   {
				   pert = rand_01() * 0.1;
				   beta *= 1. + pert;
			   }
			   accepted_moves = 0;

		   }

		   if(iter % output_log_after == 0)
		   {
			   log << "Saving best results so far " << endl;
			   save_best_results(iter);
		   }
	   }
   }



 }
//---------------------------------------------------------------------------------
void Anneal_Simulation::MPI_simulated_annealing_with_widths(int steps, double beta_start, double high_acceptance , double low_acceptance , int adapt_after, int output_log_after )
 {
   double pert;
   int accepted_moves = 0;
   double beta = beta_start;
   int accept_one = 0;
   double maxdif = 0;
   double avgdif = 0;



   if(id == 0)
   {
   	  MPI_send_pars_new();
   }
   else
   {
       MPI_receive_pars_new();
   }

   //calculate_attempt_melts(str_new);
   calculate_attempt_widths(str_new,id,proc_no);

   if(id == 0)
   {
	   MPI_receive_melts();
   }
   else
   {
	   MPI_send_melts();
   }

   if(id == 0)
   {

	   for(int i = 0; i < states.size(); i++)
	   {
		   correct_melts[i] = states[i].get_correct_width();
		   current_melts[i] = attempt_melts[i];

		   //cout << "For seq " << i << " the width is " << correct_melts[i] << endl;
		   //states[i].fill_in_predicted_yields();
	   }
	   olddifference = bestdifference = newdifference = sum_differences();

	   maxdif = get_max_difference();
	   avgdif =  get_avg_difference();

	   log << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
	   simulation_info_out << "#simbeta accratio simH AvgError MaxError  parameters " << endl;
   }

   for(int iter = 1; iter <= steps; iter++)
   {
	   accept_one = MPI_simulation_MC_step_with_widths(beta);
	   accepted_moves += accept_one;
	   if(id == 0)
	   {
		   if(accept_one)
		   {
			   maxdif = get_max_difference();
			   avgdif =  get_avg_difference();
		   }
		   log << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif<< " " << maxdif << " " << str_act << endl;

		   simulation_info_out << beta << " " <<  accepted_moves / double( (iter % adapt_after)+1) << " "  << olddifference << " " << avgdif << " " << maxdif << " " << str_act << endl;


		   if(iter % adapt_after == 0)
		   {
			   if (double(accepted_moves) / adapt_after  < low_acceptance)
			   {
				   pert = rand_01() * 0.1;
				   beta /= 1. + pert;
			   }

			   if (double(accepted_moves) / adapt_after  > high_acceptance)
			   {
				   pert = rand_01() * 0.1;
				   beta *= 1. + pert;
			   }
			   accepted_moves = 0;

		   }

		   if(iter % output_log_after == 0)
		   {
			   log << "Saving best results so far " << endl;
			   save_best_results(iter);
		   }
	   }
   }



 }
//----------------------------------------------------------------------------------
int Anneal_Simulation::MPI_simulation_MC_step_with_widths(double sim_beta)
 {
   int accepted = 0;
   if(id == 0)
   {
   	str_new.perturb_parameter(perturb_max);
   	while(!str_new.check_reasonable_limits())
   	{
   		str_new = str_act;
   		str_new.perturb_parameter(perturb_max);
   	}
   	MPI_send_pars_new();
   }
   else
   {
   	 MPI_receive_pars_new();
   }
  // cout << "Parspert " << str_new << endl;

   calculate_attempt_widths(str_new,id,proc_no);

   if(id == 0)
   {
   	MPI_receive_melts();
   }
   else
   {
   	MPI_send_melts();
   }

   if(id == 0)
   {
   	newdifference = sum_differences();
   	if ( newdifference < olddifference || rand_01()  < exp(-sim_beta * (newdifference - olddifference))   )
   	{
   		accept_move();
   		accepted = 1;

   	}
   	else
   	{
   		str_new = str_act;
   		accepted = 0;
   	}
   }
  // cout << "Parspert " << accepted << "  " << newdifference <<  " " << olddifference << endl;
   return accepted;
 }
//-------------------------------------------------------------------------------------
  int Anneal_Simulation::MPI_simulation_MC_step_with_yields(double sim_beta)
  {
    int accepted = 0;
    if(id == 0)
    {
    	str_new.perturb_parameter(perturb_max);
    	while(!str_new.check_reasonable_limits())
    	{
    		str_new = str_act;
    		str_new.perturb_parameter(perturb_max);
    	}
    	MPI_send_pars_new();
    }
    else
    {
    	 MPI_receive_pars_new();
    }
   // cout << "Parspert " << str_new << endl;

    calculate_attempt_yields(str_new,id,proc_no);

    if(id == 0)
    {
    	MPI_receive_melts();
    }
    else
    {
    	MPI_send_melts();
    }

    if(id == 0)
    {
    	newdifference = sum_differences();
    	if ( newdifference < olddifference || rand_01()  < exp(-sim_beta * (newdifference - olddifference))   )
    	{
    		accept_move();
    		accepted = 1;

    	}
    	else
    	{
    		str_new = str_act;
    		accepted = 0;
    	}
    }
   // cout << "Parspert " << accepted << "  " << newdifference <<  " " << olddifference << endl;
    return accepted;
  }
//----------------------------------------------------------------------------------------
int Anneal_Simulation::MPI_simulation_MC_step_maxdifference(double sim_beta)
{
  int accepted = 0;
  if(id == 0)
  {
	str_new.perturb_parameter(perturb_max);
	//cout << "Parspert " << str_new << endl;
	MPI_send_pars_new();
  }
  else
  {
    MPI_receive_pars_new();
  }

  calculate_attempt_melts(str_new,id,proc_no);
  if(id == 0)
  {
	  MPI_receive_melts();
  }
  else
  {
	  MPI_send_melts();
  }
  if(id == 0)
  {
	   newdifference = get_max_difference();
	   if ( newdifference < olddifference || rand_01()  < exp(-sim_beta * (newdifference - olddifference))   )
	   {
	     accept_move();
	     accepted = 1;
	   }
	   else
	   {
	     str_new = str_act;
	     accepted = 0;
	   }
	  // cout << "Parspert " << accepted << "  " << newdifference << " " << olddifference <<  endl;
  }

  return accepted;

}
//------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
int Anneal_Simulation::MPI_simulation_MC_perturb_avg_ratio_step(double sim_beta)
{
  int accepted = 0;
  if(id == 0)
  {
	str_new.perturb_avg_ratio(perturb_max);
	//cout << "Parspert " << str_new << endl;
	MPI_send_pars_new();
  }
  else
  {
    MPI_receive_pars_new();
  }

  calculate_attempt_melts(str_new,id,proc_no);
  if(id == 0)
  {
	  MPI_receive_melts();
  }
  else
  {
	  MPI_send_melts();
  }
  if(id == 0)
  {
	   newdifference = sum_differences();
	   if ( newdifference < olddifference || rand_01()  < exp(-sim_beta * (newdifference - olddifference))   )
	   {
	     accept_move();
	     accepted = 1;
	   }
	   else
	   {
	     str_new = str_act;
	     accepted = 0;
	   }
	  // cout << "Parspert " << accepted << "  " << newdifference << " " << olddifference <<  endl;
  }

  return accepted;

}
//-----------------------------------------------------------------------------------------
int Anneal_Simulation::MPI_simulation_STUN_MC_step(double sim_beta)
{
  int accepted = 0;
  if(id == 0)
  {
	  //str_new.perturb_parameter(perturb_max);
	  int notall_in_limits = 0;
	  do
	  {
		  str_new.perturb_parameter(perturb_max);
		  notall_in_limits = 0;
		  for(int i = 0; i < 2; i++)
		  {
			  if(str_new.vals[i] < 0.7 || str_new.vals[i] > 2.2)
			  {
				  notall_in_limits = 1;
				  //break;
			  }
		  }
		  for(int i =2; i < PAR_SIZE; i++)
		  {
			  if(str_new.vals[i] < 0.8 || str_new.vals[i] > 3.5)
			  {
				  notall_in_limits = 1;
				  //break;
			  }
		  }
		  if(notall_in_limits)
			  str_new = str_act;

	  } while(notall_in_limits);
	//cout << "Parspert " << str_new << endl;
	MPI_send_pars_new();
  }
  else
  {
    MPI_receive_pars_new();
    //cout << "Received " << str_new << endl;
  }

  calculate_attempt_melts(str_new,id,proc_no);
  if(id == 0)
  {
	  MPI_receive_melts();
  }
  else
  {
	  MPI_send_melts();
  }
  if(id == 0)
  {
	     accepted = 0;

	    // cout << "Parspert " << str_new << endl;
	     double _gama = 10.;
	     //calculate_attempt_melts(str_new);
	     newdifference = sum_differences();

	     double e_stun_new = 1.0 - exp(-_gama * (newdifference - bestdifference) );
	     double e_stun_old = 1.0 - exp(-_gama * (olddifference - bestdifference) );

	     if ( newdifference < olddifference || rand_01()  < exp(-sim_beta * (e_stun_new - e_stun_old))   )
	     {
	       accept_move();
	       accepted = 1;

	     }
	     else
	     {
	       str_new = str_act;
	       accepted = 0;
	     }

	  // cout << "Parspert " << accepted << "  " << newdifference << " " << olddifference <<  endl;
  }

  return accepted;

}
//------------------------------------------------------------------------------------------
int Anneal_Simulation::MPI_receive_melts(void)
{
 if(id != 0)
	 return 0;

 int TAG = 0;
 MPI_Status stat;
 int seq_per_proc = states.size() / proc_no;
 int stop ;
 int bufsize = seq_per_proc;
 int ret = 0;
 for(int node_id = 1; node_id < proc_no; node_id++)
 {
	 bufsize = seq_per_proc;
	 stop =  (node_id + 1) *  (states.size() / proc_no );
	 if(node_id == proc_no - 1)
	 {
	 	  stop = states.size();
	 	  bufsize = stop -  node_id * (seq_per_proc) ;
	 }
     ret = MPI_Recv( (void *) (attempt_melts+node_id*seq_per_proc),bufsize,MPI_DOUBLE,node_id,TAG,MPI_COMM_WORLD,&stat);
     if(ret != MPI_SUCCESS)
     {
    	 throw string("Error, did not receive correctly melting temperatures");
     }
 }
 return 1;
}

//----------------------------------------------------------------------------------------
int Anneal_Simulation::MPI_send_melts(void)
{
  int stop =  (id + 1) *  (states.size() / proc_no );
  int seq_per_proc = states.size() / proc_no;
  int bufsize = seq_per_proc;
  int TAG = 0;
  if(id == proc_no - 1)
  {
	  stop = states.size();
	  bufsize = stop -  id * (seq_per_proc) ;
  }
  int ret;
  if(id != 0)
  {
	ret =  MPI_Send((void *)( attempt_melts+id*seq_per_proc),bufsize,MPI_DOUBLE,0,TAG,MPI_COMM_WORLD);
    /*
	for(int k = 0; k < bufsize; k++)
	{
		cout << "Send " << beta_to_kelvin(attempt_melts[id*seq_per_proc+k]) << " " << endl;
	}
	*/
  }
  else return 0;

  if(ret != MPI_SUCCESS)
  {
	  throw string("Not successfully sent melting temeperatures");
  }

  return 1;
}

//---------------------------------------------------------------------------------------------
int Anneal_Simulation::MPI_send_pars_new(void)
{
   int TAG = 1;
   int bufsize = PAR_SIZE;
   int ret;
   if(id == 0) //I am sending to all the nodes
   {
     for(int node_id = 1; node_id < proc_no; node_id++)
     {
      ret = MPI_Send((void *)str_new.vals,bufsize,MPI_DOUBLE,node_id,TAG,MPI_COMM_WORLD);
      //cout << "Sending from " << id << " " << str_new << endl;
      if(ret != MPI_SUCCESS)
     	    	   throw string("Error while receiving new values");
     }
   }
   else return 0;

   return 1;
}


//--------------------------------------------------------------------------------------
int Anneal_Simulation::MPI_receive_pars_new(void)
{
	 int TAG = 1;
	 int bufsize = PAR_SIZE;
	 MPI_Status stat;
     int ret;
     if(id != 0)
	 {

	      ret =  MPI_Recv((void *)str_new.vals,bufsize,MPI_DOUBLE,0,TAG,MPI_COMM_WORLD,&stat);
	      //cout << "Receiving in " << id <<" " <<  str_new << endl;
	      if(ret != MPI_SUCCESS)
	    	   throw string("Error while receiving new values");

	 }
	 else return 0;

  return 1;
}
//-----------------------------------------------------------------------------------
#endif
