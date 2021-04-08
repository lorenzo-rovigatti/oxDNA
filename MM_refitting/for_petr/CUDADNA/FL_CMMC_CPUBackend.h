/*
 * FL_CMMC_CPUBackend.h
 *
 *  Created on: 26/nov/2010
 *      Author: flavio 
 */

#ifndef FL_CMMC_CPUBACKEND_H_
#define FL_CMMC_CPUBACKEND_H_

#include "MCBackend.h"
#include "MC_CPUBackend.h"
#include "Weights.h"
#include "OrderParameters.h"
#include "Histogram.h"

template<typename number>
class FL_CMMC_CPUBackend: public MC_CPUBackend<number> {
protected:
	Particle<number> *_particles_old;

	bool _have_us;
	bool _reload_hist;
	OrderParameters _op;
	char _op_file[512];
	Weights _w;
	char _weights_file[512];
	Histogram _h;
	char _traj_hist_file[512];
	char _last_hist_file[512];
	char _init_hist_file[512];
	char _state_str[512];

	int _maxclust;
	number _e_min_cutoff;
	number _e_max_cutoff;

	number _e_stack;

	bool _adjust_delta;
	double _my_delta_tr;
	double _my_delta_or;
	llint _acc_p_tr, _acc_p_or, _eq_steps, _adj_interval;
	llint _try_or_p, _try_tr_p;

	// stuff for Frenkel-Ladd
	bool _have_fl;
	number _eta, _lambda_or, _lambda_tr;
	int _ref_index;
	Particle<number> * _ref_particles;
	LR_vector<number> _cdm;
	char _fl_str[512];
	
	int _netemps;
	double * _etemps;
	
	inline number _excluded_volume(const LR_vector<number> &r, number sigma, number rstar, number b, number rc);

	inline number _particle_particle_bonded_interaction_n5 (Particle<number> *p, Particle<number> *q);
	inline number _particle_particle_bonded_interaction_n3 (Particle<number> *p, Particle<number> *q);
	inline number _particle_particle_interaction_pq (Particle<number> *p, Particle<number> *q);
	inline number _particle_particle_interaction (Particle<number> *p, Particle<number> *q);
	inline number _particle_particle_hb (Particle<number> *p, Particle<number> *q);
	inline number _compute_stacking_energy ();

	// functions for FL
	inline number _particle_einst_tr_energy (Particle<number> *p);
	inline number _particle_einst_tr_energy_old (Particle<number> *p);
	inline number _particle_einst_tr_delta_energy (Particle<number> *p);
	inline number _particle_einst_tr_delta_energy_old (Particle<number> *p);
	inline number _cluster_einst_tr_delta_energy (int *, int);
	inline number _cluster_einst_tr_delta_energy_old (int *, int);
	inline number _particle_einst_or_energy (Particle<number> *p);
	inline number _particle_einst_energy (Particle<number> *p);
	inline number _system_einst_energy ();
	inline number _system_einst_energy_old ();
	inline number _system_einst_tr_energy ();
	inline number _system_einst_tr_energy_old ();
	inline number _system_einst_or_energy ();
	
	inline void build_cluster (int seed, number cutoff, int maxsize, int * clust, int * size);
	inline bool check_cluster (int seed, number cutoff, int * clust, int nclust, number * ecluster);
	inline void _translate_cluster (int *, int);
	inline void _rotate_particle(Particle<number> *p);
	inline void _rotate_cluster (int *, int);
	inline number _cluster_energy(int *, int, bool);
	inline void _fill_e_neigh (Particle<number> *, Particle<number> *, number, int);
	inline void _fill_h_bonds (Particle<number> *, Particle<number> *, bool);
	void _update_metainfo ();
	void _update_ops ();

public:
	FL_CMMC_CPUBackend(IOManager *IO);
	virtual ~FL_CMMC_CPUBackend();

	virtual void get_settings (input_file &inp);
	virtual void print_conf (llint curr_step, bool reduced, bool only_last);
	virtual void print_conf (llint curr_step, bool only_last);
	virtual void print_energy (llint curr_step);
	
	number get_e_stack () const { return _e_stack; }
	
	//void init(ifstream &conf_input);
	void init(char conf_filename[256]);
	
	void init_fl();

	void sim_step(llint cur_step);
	inline void check_energy();
	inline void check_ops();
	char * get_fl_str ();
	char * get_op_state_str ();
};

#endif /* FL_CMMC_CPUBACKEND_H_ */
