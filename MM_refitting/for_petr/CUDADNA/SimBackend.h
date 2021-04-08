/*
 * SimBackend.h
 *
 *  Created on: 03/set/2010
 *      Author: lorenzo
 */

#ifndef SIMBACKEND_H_
#define SIMBACKEND_H_

#define SIM_MD 0
#define SIM_MC 1

#include <cmath>
#include <fstream>
#include <float.h>

extern "C" {
#include "parse_input/parse_input.h"
#include "timing/timing.h"
}

#include "defs.h"
#include "Utils.h"
#include "Particle.h"
#include "Interaction.h"

using namespace std;

class IOManager;

class ISimBackend {
public:
	ISimBackend() {};
	virtual ~ISimBackend() {};
	
	llint _start_step_from_file;

	virtual void get_settings(input_file &inp) = 0;
	//virtual void init(ifstream &conf_input) = 0;
	virtual void init(char path[256]) = 0;
	virtual void sim_step(llint curr_step) = 0;
	virtual void print_energy(llint curr_step) = 0;
	virtual void print_conf(llint curr_step, bool reduced=false, bool only_last=false) = 0;
	virtual int get_N_updates() = 0;
};

template<typename number>
class SimBackend: public ISimBackend{
	friend class IOManager;
protected:
	// this pointer makes possible to use the forward declaration
	IOManager *_IO;

	LR_timer _timer;
	int _timer_msgs_number;
	msg _timer_msgs[10];

	bool _print_timings;
	char _timings_filename[256];
	bool _is_CUDA_sim;

	char _topology_filename[256];

	bool _external_forces;
	char _external_filename[256];

	bool _grooving;

	int _sim_type;

	int _N;
	int _N_strands;
	int _max_neigh;
	number _T;
	number _box_side;
	number _verlet_skin;
	number _sqr_verlet_skin;
	number _sqr_rverlet;

	// interaction's parameters
	Interaction<number> _interaction;
	number _rcut;
	number _sqr_rcut;

	// potential and kinetic energies
	number _U, _K;
	
	// potential energy due to hydrogen bonding
	number _U_hydr;
	
	// potential energy due to hydrogen bonding
	number _U_stack, _dU_stack;
       
        //change in potential energy
	number _dU; 

	Particle<number> *_particles;

	void _get_number_settings(input_file &inp);

	// cells and lists management
	bool _are_lists_old;
	int *_heads;
	int _N_cells;
	int _N_cells_side;
	int _N_updates;
	int _confs_to_skip;

	void _check_input_sanity();
	void _read_topology();
	void _read_external_forces ();
	void _fill_cell_index(const LR_vector<number> &pos, int cell_index[3]);
	void _create_cells();
	void _update_lists();

public:
	SimBackend(IOManager *IO);
	virtual ~SimBackend();

	virtual void get_settings(input_file &inp);
	//virtual void init(ifstream &conf_input);
	virtual void init(char path[256]);

	int get_N_updates() {return _N_updates; }
};

#endif /* SIMBACKEND_H_ */
