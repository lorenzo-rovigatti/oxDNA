/**
 * @file    MD_MPIBackend.h
 * @date    31/12/2011
 * @author  petr
 *
 *
 */

#ifndef MD_MPIBACKEND_H_
#define MD_MPIBACKEND_H_

#include "MD_CPUBackend.h"

#include <mpi.h>
//#include <mpi>


struct Serialized_particle_position
{
	int index;
	LR_vector pos;
	/// these are the positions of backbone, base and stack sites relative to the center of mass
	LR_vector pos_back;
	LR_vector pos_stack;
	LR_vector pos_base;
	LR_matrix orientation;
	/// transpose (= inverse) orientational matrix
	LR_matrix orientationT;
	int _N_neigh;
	void read_from(BaseParticle &par);
	int write_to(BaseParticle &par);
};

struct Serialized_particle_force_torque
{
	int index;
	LR_vector torque;
	LR_vector force;

	Serialized_particle_force_torque(int _index = -1)
	{torque = LR_vector(0,0,0); force = LR_vector(0,0,0); index = _index;}

	void read_from(BaseParticle &par);
	int add_to(BaseParticle &par);

};

struct Master_to_node_info {

	int min_index;
	int max_index;
	int particle_count;

	Master_to_node_info(int _min=0, int _max=0, int _count=0) :  min_index(_min), max_index(_max), particle_count(_count) {}
    void set_no_change(void) {particle_count = -1;}
    int no_change(void) {return particle_count == -1 ? 1 : 0; }
};

 struct Energy_info
{
 number U_hydr;
 number U;
 Energy_info(number _U = 0, number _U_hydr =0 ) {U = _U; U_hydr = _U_hydr;}
};




class MD_MPIBackend: public MD_CPUBackend {
protected:
	int _myid; ///id of the machine
	int _proc_size; ///total number of cores in the simulation
    int _number_of_interactions; ///number of interaction to be treated by this processor
   // int *_particles_to_process; //list of particle indices that will be processed by this processor
   // int *_communicate_particles; //communicate_particles is a list of particles to process, which is sent by master node

    Serialized_particle_position *_serialized_particles;
    Serialized_particle_force_torque *_serialized_forces;
    Master_to_node_info _info_process;
    int _total_no_of_interactions;
    int _interaction_per_node;
    bool _mpi_lists_are_old;


    int _MPI_send_block_data(void *data,size_t size,int node_to,int TAG=1);
    int _MPI_receive_block_data(void *data, size_t size, int node_from, int TAG=1);

	//number _MPI_particle_particle_bonded_interaction(BaseParticle *p);
	//number _MPI_particle_particle_interaction(BaseParticle *p, BaseParticle *q);

    void _Serialize_all_particles(void);
    void _Serialize_all_forces(void);

    int _MPI_send_master_to_slave_info(Master_to_node_info& inf,int send_to);
    int _MPI_receive_master_to_slave_info(Master_to_node_info& inf,int rec_from=0);

	//int _MPI_send_all_serialized_particles_to_slave(int slave_id); //sends a particle with all its neighbors
	int _MPI_send_serialized_particle_to_slave(Serialized_particle_position& part,int slave_id);
    int _MPI_send_serialized_particle_with_neighbors_to_slave(Serialized_particle_position &part,int slave_id);
    int _MPI_send_all_serialized_particles_with_neighbors_to_slave(int slave_id);


    int _MPI_receive_serialized_particle(int from_id=0);
    int _MPI_receive_serialized_particle_with_neighbors(int from_id = 0);

    int _MPI_receive_and_fill_forces_from_slave(int slave_id);
    int _MPI_send_serialized_forces_to_master(int master_id = 0);



	int _MPI_master_send_interactions(void);
	int _MPI_slave_receive_interactions(void);

	void _Evaluate_my_particles(void); //each node evaluates only the particles in this reposnsible for
	void _MPI_compute_forces();
	//void _MPI_update_lists();

public:
	MD_MPIBackend();
	virtual ~MD_MPIBackend();

	virtual void init();
    virtual void print_conf(bool reduced=false, bool only_last=false) ;
	virtual void sim_step();

	//void activate_thermostat();
};

#endif /* MD_MPIBACKEND_H_ */
