/*
 * MD_MPIBackend.cpp
 *
 *  Created on: 1/1/2012
 *      Author: petr
 */

#include "MD_MPIBackend.h"


//---------------------------------------------------------------------------------------
 void Serialized_particle_force_torque::read_from(BaseParticle &par)
{
	index = par.index;
	torque = par.torque;
	force = par.force;
}
//----------------------------------------------------------------------------------------
 int Serialized_particle_force_torque::add_to(BaseParticle &par)
{
	if(par.index != index)
	{
		cerr << "CRITICAL error, index incompatibility, attempting to overwrite force particle id " << par.index << " with " << index << endl;
		MPI_Abort(MPI_COMM_WORLD,1);
		return -1;
	}
	par.torque += torque;
	par.force += force;

	return index;
}
//---------------------------------------------------------------------------------------
 void Serialized_particle_position::read_from(BaseParticle &par)
{
	index = par.index;
	pos = par.pos;
	pos_back = par.pos_back;
	pos_stack = par.pos_stack;
	pos_base = par.pos_base;
	orientation= par.orientation;
	orientationT = par.orientationT;
	_N_neigh = par._N_neigh;
}
//----------------------------------------------------------------------------------------
 int Serialized_particle_position::write_to(BaseParticle &par)
{
	if(par.index != index)
	{
		cerr << "CRITICAL error, index incompatibility, attempting to overwrite position particle id " << par.index << " with " << index << endl;
		MPI_Abort(MPI_COMM_WORLD,1);
		return -1;

	}

	par.pos = pos;
	par.pos_back = pos_back;
	par.pos_stack = pos_stack;
	par.pos_base = pos_base;
	par.orientation= orientation;
	par.orientationT = orientationT;
	par._N_neigh = _N_neigh;
	return par.index;
}
//---------------------------------------------------------------------------------------

int MD_MPIBackend::_MPI_send_block_data(void *data,size_t size,int node_to,int TAG)
{
	  int ret =  MPI_Send((void *)data,size,MPI_CHAR,node_to,TAG,MPI_COMM_WORLD);
	  if(ret != MPI_SUCCESS)
	  {
		   	throw oxDNAException("Error while sending MPI message");
	  }
	  else return 1;

	  return 0;
}
//---------------------------------------------------------------------------------------

int MD_MPIBackend::_MPI_receive_block_data(void *data, size_t size, int node_from, int TAG)
{
	     MPI_Status stat;
	     int ret = MPI_Recv( (void *) data,size,MPI_CHAR,node_from,TAG,MPI_COMM_WORLD,&stat);
	     if(ret != MPI_SUCCESS)
	     {
	    	throw oxDNAException("Error while receving MPI message");

	     }
	     else return 1;

	     return 0;
}

//---------------------------------------------------------------------------------------

MD_MPIBackend::MD_MPIBackend() : MD_CPUBackend() {
	this->_mpi_lists_are_old = true;

	MPI_Comm_rank (MPI_COMM_WORLD, &(this->_myid));
	MPI_Comm_size (MPI_COMM_WORLD, &(this->_proc_size));

	//cout << "Initialized process " << this->_myid <<  " of " << this->_proc_size << endl;
//	this->_particles_to_process = new int[this->_N+1];
//	this->_communicate_particles = new int[this->_N+1];

	/*
	MPI_Datatype types[] = { MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_INTEGER, MPI_FLOAT };
	std::vector<int> len(6, 1);
	std::vector<MPI_Aint> disp(6, 0);
	particle temp;
	MPI_Aint base;
	MPI_Address(&temp, &base);
	MPI_Address(&temp._x, &disp[0]);
	MPI_Address(&temp._y, &disp[1]);
	MPI_Address(&temp._xvel, &disp[2]);
	MPI_Address(&temp._yvel, &disp[3]);
	MPI_Address(&temp._isStaticInt, &disp[4]);
	MPI_Address(&temp._size, &disp[5]);
	for (int i=0; i<6; ++i)
	{
		disp[i] = disp[i] - base;
	}
	MPI_Type_struct(6, &len[0], &disp[0], types, &_particleType);
	MPI_Type_commit(&_particleType);
*/
}

//-----------------------------------------------------------------------------------------


MD_MPIBackend::~MD_MPIBackend() {

	//delete [] this->_particles_to_process;
	//delete [] this->_communicate_particles;
	delete [] this->_serialized_particles ;
	delete [] this->_serialized_forces ;
}
//-------------------------------------------------------------------------------------------------------
int MD_MPIBackend::_MPI_send_serialized_particle_to_slave(Serialized_particle_position& part,int slave_id)
{
	return this->_MPI_send_block_data( (void *)&part,sizeof(part),slave_id);
}
//---------------------------------------------------------------------------------------------------

int MD_MPIBackend::_MPI_send_serialized_particle_with_neighbors_to_slave(Serialized_particle_position &part,int slave_id)
{
	 this->_MPI_send_block_data( (void *)&part,sizeof(part),slave_id);
 	 int ret =  this->_MPI_send_block_data( static_cast<void *>(this->_particles[part.index].get_verlet_list()) ,sizeof(int)*part._N_neigh,slave_id);
     return ret;
}
//---------------------------------------------------------------------------------------------------
int MD_MPIBackend::_MPI_receive_serialized_particle(int from_id)
{
	Serialized_particle_position part;
	this->_MPI_receive_block_data( (void *)&part,sizeof(part),from_id);
	if(part.index < 0 || part.index >= this->_N)
		throw oxDNAException(" Received particle with impossible index %d",part.index);
    part.write_to(this->_particles[part.index]);
	return part.index;
}
//---------------------------------------------------------------------------------------------------
int MD_MPIBackend::_MPI_receive_serialized_particle_with_neighbors(int from_id)
{
	Serialized_particle_position part;
	this->_MPI_receive_block_data( (void *)&part,sizeof(part),from_id);

	if(part.index < 0 || part.index >= this->_N)
		throw oxDNAException(" Received particle with impossible index %d",part.index);

	this->_MPI_receive_block_data( static_cast<void *>(this->_particles[part.index].get_verlet_list()),sizeof(int)*part._N_neigh,from_id);
	part.write_to(this->_particles[part.index]);
	return part.index;
}
/*
//----------------------------------------------------------------------------------------------------
int MD_MPIBackend::_MPI_receive_force_from_slave(int slave_id)
{
	Serialized_particle_force_torque ftrq;
	_MPI_receive_block_data( (void *)&ftrq,sizeof(ftrq),slave_id);

	if(ftrq.index < 0 || ftrq.index >= this->_N)
			throw oxDNAException(" Received particle with impossible index %d",part.index);

	ftrq.add_to(this->_particles[ftrq.index]);
	return ftrq.index;
}
//---------------------------------------------------------------------------------------------------
int MD_MPIBackend::_MPI_send_force_to_master(int particle_id,int master_id)
{
	Serialized_particle_force_torque ftrq;
	ftrq.read_from(this->_particles[particle_id]);
	return _MPI_send_block_data( (void *)&ftrq,sizeof(ftrq),master_id);
}
*/
//---------------------------------------------------------------------------------------------------
int MD_MPIBackend::_MPI_send_master_to_slave_info(Master_to_node_info& inf,int send_to)
{
   return this->_MPI_send_block_data( (void *)&inf,sizeof(inf),send_to);
}
//-----------------------------------------------------------------------------------------------------
int MD_MPIBackend::_MPI_receive_master_to_slave_info(Master_to_node_info& inf,int rec_from)
{
   return this->_MPI_receive_block_data( (void *)&inf,sizeof(inf),rec_from);
}

//--------------------------------------------------------------------------------------------------------
void MD_MPIBackend::_Serialize_all_particles(void)
{
  for(int i = 0; i < this->_N; i++)
  {
    this->_serialized_particles[i].read_from(this->_particles[i]);
  }
}
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MD_MPIBackend::_Serialize_all_forces(void)
{
  for(int i = 0; i < this->_N; i++)
  {
    this->_serialized_forces[i].read_from(this->_particles[i]);
  }
}
//----------------------------------------------------------------------------------------------------------
void MD_MPIBackend::_Evaluate_my_particles(void)
{

	int neigh;
	BaseParticle *p;

	//forces have to be set to 0 (on slave nodes) before evaluating this:
	this->_U = this->_U_hydr = (number) 0;
	for(int i = this->_info_process.min_index; i <= this->_info_process.max_index; i++) {
		p = &this->_particles[i];
		this->_U += _particle_particle_bonded_interaction(p);

		std::vector<BaseParticle *> neighs = this->_lists->get_neigh_list(p);
		for(unsigned int n = 0; n < neighs.size(); n++) {
			BaseParticle *q = neighs[n];
			this->_U += this->_interaction->pair_interaction_nonbonded(p, q, true, true);
		}
	}
	//this can be optimized
	if(this->_myid != 0)
	{
		this->_Serialize_all_forces();
	}

}
//---------------------------------------------------------------------------------------------------------


void MD_MPIBackend::sim_step() {

  if(this->_myid == 0)   //this is master node
  {
    get_time(&this->_timer, 0);

	get_time(&this->_timer, 2);
	MD_CPUBackend::_first_step();
	get_time(&this->_timer, 3);

	get_time(&this->_timer, 6);
	if(!this->_lists->is_updated()) {
		this->_lists->global_update();
		this->_N_updates++;
		this->_mpi_lists_are_old = true;
	}
	get_time(&this->_timer, 7);

	get_time(&this->_timer, 8);

	this->_MPI_compute_forces();

	MD_CPUBackend::_second_step();

	get_time(&this->_timer, 9);

	get_time(&this->_timer, 10);
	_thermostat->apply(_particles, current_step());
	get_time(&this->_timer, 11);

	get_time(&this->_timer, 1);

	process_times(&this->_timer);
  }
  else
  {
  	this->_MPI_compute_forces();  //slave node just computes forces
  }
}
//------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------
int MD_MPIBackend::_MPI_send_all_serialized_particles_with_neighbors_to_slave(int slave_id)
{
  if(slave_id > this->_proc_size - 1)
  {
	  throw oxDNAException("Trying to send data to nonexistent node %d",slave_id);

  }

  int retval=1;

  for(int i = 0; i < this->_N; i++)
  {
	retval = this->_MPI_send_serialized_particle_with_neighbors_to_slave(this->_serialized_particles[i],slave_id);
  }

  return retval;
}
//-------------------------------------------------------------------------------------------
int MD_MPIBackend::_MPI_master_send_interactions(void) {
	Master_to_node_info info_struct_send;
	if(this->_mpi_lists_are_old) //we are refreshing neighbour lists
	{
	    this->_total_no_of_interactions = 0;

		for(int i = 0; i < this->_N; i++)
	    {
			if(this->_particles[i].n5 != P_VIRTUAL)
				this->_total_no_of_interactions++;
			this->_total_no_of_interactions += this->_particles[i].get_N_neigh();
			this->_serialized_particles[i].read_from(this->_particles[i]); //serialization of data
	    }

		this->_interaction_per_node = (this->_total_no_of_interactions / this->_proc_size) + 1;
        //cout << "Int per node is " << this->_interaction_per_node << " out of " << this->_total_no_of_interactions << endl;
		int interactions_filled = 0;
		int min_index = 0;
		int max_index = 0;
		int process_to_send_data = 0;
		int interactions_per_node_filled = 0;
		for(int i = 0; i < this->_N; i++)
		{

		  interactions_filled += this->_particles[i].get_N_neigh();
		  interactions_per_node_filled  += this->_particles[i].get_N_neigh();
		  if(this->_particles[i].n5 != P_VIRTUAL)
		  {
			  interactions_filled++;
			  interactions_per_node_filled++;
		  }
		  if(interactions_per_node_filled >= this->_interaction_per_node || i == this->_N-1)
		  {
			max_index = i;
			process_to_send_data++;
			if(process_to_send_data >= this->_proc_size) //this is just for evaluation by the master node 0
			{
				//by now, all interactions have been filled
				if(interactions_filled != this->_total_no_of_interactions || max_index != this->_N-1)
				{
					throw oxDNAException(" Error in MPI data distribution, bug in communication master->slave enumeration of states, filled %d out of %d ",interactions_filled,this->_total_no_of_interactions );
				}
				//cout << "Process 0 evaluates from " << min_index << " to " << max_index << endl;
				//this->_IO->log(0,"Process 0 evaluates from %d to %d",min_index,max_index);

				this->_info_process.min_index = min_index;  //nodes that master is reponsible fot
				this->_info_process.max_index = max_index;
			}
			else  //in this case, evaluation will be done by slave node
			{
				//cout << "Process " << process_to_send_data << "  evaluates from " << min_index << " to " << max_index << endl;
				info_struct_send.min_index = min_index;
				info_struct_send.max_index = max_index;
				info_struct_send.particle_count = this->_N; //in this case everything is sent
				this->_MPI_send_master_to_slave_info(info_struct_send,process_to_send_data);
				for(int part = 0; part < this->_N; part++) //sending everything
				{
				  //cout << " Sending " << this->_serialized_particles[part].index << endl;
				  this->_MPI_send_serialized_particle_with_neighbors_to_slave(this->_serialized_particles[part],process_to_send_data);
				  //cout << " Sent serialized particles with neighbors " << endl;
				}
			}
			interactions_per_node_filled = 0;
			min_index = max_index+1;

		  }
		}

		if(process_to_send_data < this->_proc_size  ) //we did not fill all processes with work
		{
			process_to_send_data++;
			info_struct_send.min_index = 1;
			info_struct_send.max_index = 0;
			info_struct_send.particle_count = 0;
			while(process_to_send_data < this->_proc_size - 1)
			{
			  this->_MPI_send_master_to_slave_info(info_struct_send,process_to_send_data);
			  process_to_send_data++;
			}
			this->_info_process.min_index = 1;
			this->_info_process.max_index = 0;
			this->_info_process.particle_count = 0;

		}

		this->_mpi_lists_are_old = false;
	}
	else //no list refreshing, lets just send all
	{
	  this->_Serialize_all_particles();
      info_struct_send.set_no_change();
      for(int proc_id = 1; proc_id < this->_proc_size; proc_id++)
      {
    	  this->_MPI_send_master_to_slave_info(info_struct_send,proc_id);
    	  for(int part = 0; part < this->_N; part++) //sending everything
    	  {
    	    this->_MPI_send_serialized_particle_to_slave(this->_serialized_particles[part],proc_id);
    	    //cout << "Sent serialized particles without neighbors" << endl;
    	  }

      }
	}

 return 1;
}
//--------------------------------------------------------------------------------------------------
int MD_MPIBackend::_MPI_slave_receive_interactions(void)
{
	Master_to_node_info info;

	this->_MPI_receive_master_to_slave_info(info);
	if(info.no_change()) //no neighbor list change
	{
		for(int i =0; i < this->_N; i++)
		{
			  this->_MPI_receive_serialized_particle();
		}
		//cout << "Received particles without neighbors" << endl;
	}
	else
	{
      this->_info_process = info;
      if(this->_info_process.particle_count > 0)
      {
    	  for(int i =0; i < this->_N; i++)
    	  {
    		   this->_MPI_receive_serialized_particle_with_neighbors();
    		   //cout << "Received particles with neighbors" << endl;
    	  }
      }
	}

	//now, we set all forces in existing particles to 0
	for(int i =0; i < this->_N; i++)
	{
	  this->_particles[i].force = LR_vector(0,0,0);
	  this->_particles[i].torque = LR_vector(0,0,0);
	}

	return info.particle_count;
}
//--------------------------------------------------------------------------------------------------

int MD_MPIBackend::_MPI_receive_and_fill_forces_from_slave(int slave_id)
{
	Energy_info einfo;
	_MPI_receive_block_data((void *)(this->_serialized_forces), sizeof( Serialized_particle_force_torque)*this->_N, slave_id);
	int retval = _MPI_receive_block_data((void *)&einfo,sizeof(einfo),slave_id);

	this->_U += einfo.U;
	this->_U_hydr += einfo.U_hydr;

	for(int i = 0; i < this->_N; i++)
	{
		this->_serialized_forces[i].add_to(this->_particles[i]);
	}

	return retval;
}
//--------------------------------------------------------------------------------------------------

int MD_MPIBackend::_MPI_send_serialized_forces_to_master(int master_id)
{
	 Energy_info einfo;
	 einfo.U = this->_U;
	 einfo.U_hydr = this->_U_hydr;
	 this->_MPI_send_block_data((void *)(this->_serialized_forces), sizeof( Serialized_particle_force_torque)*this->_N, master_id);
	 int retval = this->_MPI_send_block_data((void *)&einfo,sizeof(einfo),master_id);

	 return retval;
}

//---------------------------------------------------------------------------------------------------


void MD_MPIBackend::_MPI_compute_forces() {
 if(this->_myid == 0)
 {
	 this->_MPI_master_send_interactions();
	 this->_Evaluate_my_particles();
	 for(int sl_id = 1; sl_id < this->_proc_size; sl_id++)
	 {
	   this->_MPI_receive_and_fill_forces_from_slave(sl_id);
	 }
 }
 else //slave node
 {
	 this->_MPI_slave_receive_interactions();
	 this->_Evaluate_my_particles();
	 this->_MPI_send_serialized_forces_to_master();
 }

}
/*
//-------------------------------------------------------------------

void MDBackend::print_conf(llint curr_step, bool reduced, bool only_last) {
	if(this->_myid == 0) {
	 if(reduced == false) this->_IO->print_conf(*this, curr_step, only_last);
	 else this->_IO->print_reduced_conf(*this, curr_step);
	}
}
*/
//-----------------------------------------------------------------

//void MD_MPIBackend::init(ifstream &conf_input) {
void MD_MPIBackend::init() {
//	MD_CPUBackend::init(conf_input);
	MD_CPUBackend::init();
	this->_serialized_particles = new Serialized_particle_position[this->_N];
	this->_serialized_forces = new Serialized_particle_force_torque[this->_N];
	//cout << "INIT finished, this->N is " << this->_N << endl;

}
//---------------------------------------------------------------------------

void MD_MPIBackend::print_conf(bool reduced, bool only_last) {
	if(this->_myid == 0) SimBackend::print_conf(reduced, only_last);
}
//--------------------------------------------------------------------------
//----------------------------------------------------------------------
