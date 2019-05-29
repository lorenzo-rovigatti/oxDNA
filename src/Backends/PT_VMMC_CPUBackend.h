/**
 * @file    PT_VMMC_CPUBackend.h
 * @date    26/sep/2012
 * @author  flavio
 *
 *
 */

#ifndef _PT_VMMC_CPUBACKEND_H_
#define _PT_VMMC_CPUBACKEND_H_

#include "VMMC_CPUBackend.h"
#include "mpi.h"

template<typename number>
struct PT_serialized_particle_info {
	LR_vector<number> pos;
	LR_matrix<number> orientation;
	// transpose (= inverse) orientational matrix
	void read_from(BaseParticle<number> * par);
	void write_to(BaseParticle<number> * par);
};

template <typename number>
struct PT_energy_info {
	number U, U_hydr, U_stack, U_ext, T, w;
	int weight_index, replica_id;
	PT_energy_info (number _U = 0, number _U_hydr = 0, number _U_stack = 0., number _T = 1.) {
		U = _U;
		U_hydr = _U_hydr;
		U_stack = _U_stack;
		T = _T;
		w = (number) 1.;
		weight_index = 0;
		U_ext = (number) 0.;
		replica_id = 0;
	}
};

template<typename number>
class PT_VMMC_CPUBackend: public VMMC_CPUBackend<number> {
protected:
	int _npttemps;
	double * _pttemps;
	bool _pt_common_weights;

	char _replica_info[256];

	/// keep track of who's got which replica
	int _which_replica;

	number _U_ext;

	llint _pt_exchange_tries, _pt_exchange_accepted;

	Weights _irresp_w;
	char _irresp_weights_file[512];

	int _pt_move_every;
	int _my_mpi_id, _mpi_nprocs;
	int _MPI_send_block_data(void *data, size_t size, int node_to,int TAG=1);
	int _MPI_receive_block_data(void *data, size_t size, int node_from, int TAG=1);

	void _build_exchange_conf();
	void _rebuild_exchange_conf();
	void _build_exchange_energy();
	void _rebuild_exchange_energy();
	void _get_exchange_conf (int other);
	void _get_exchange_energy (int other);
	void _send_exchange_conf (int other);
	void _send_exchange_energy (int other);

	PT_serialized_particle_info<number> * _exchange_conf;
	PT_energy_info<number> _exchange_energy;

	bool _oxDNA2_stacking;

public:
	PT_VMMC_CPUBackend();
	virtual ~PT_VMMC_CPUBackend();
	//void init(ifstream &conf_input);
	void init();
	void get_settings (input_file &inp);

	int get_mpi_id () { return _my_mpi_id; }
	int get_which_replica () { return _which_replica; }
	char * get_replica_info_str ();
	number get_pt_acc ();

	void sim_step(llint cur_step);
};

#endif /* PT_VMMC_CPUBACKEND_H_ */

