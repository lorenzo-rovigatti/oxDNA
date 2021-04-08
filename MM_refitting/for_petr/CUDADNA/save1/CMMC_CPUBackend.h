/*
 * CMMC_CPUBackend.h
 *
 *  Created on: 26/nov/2010
 *      Author: flavio 
 */

#ifndef CMMC_CPUBACKEND_H_
#define CMMC_CPUBACKEND_H_

#include "MCBackend.h"
#include "MC_CPUBackend.h"
#include "Weights.h"
#include "OrderParameters.h"
#include "Histogram.h"

template<typename number>
class CMMC_CPUBackend: public MC_CPUBackend<number> {
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
	
	int _netemps;
	double * _etemps;
	
	inline number _excluded_volume(const LR_vector<number> &r, number sigma, number rstar, number b, number rc);

	inline number _particle_particle_bonded_interaction_n5 (Particle<number> *p, Particle<number> *q);
	inline number _particle_particle_bonded_interaction_n3 (Particle<number> *p, Particle<number> *q);
	inline number _particle_particle_interaction_pq (Particle<number> *p, Particle<number> *q);
	inline number _particle_particle_interaction (Particle<number> *p, Particle<number> *q);
	inline number _particle_particle_hb (Particle<number> *p, Particle<number> *q);
	inline number _stacking_energy_n3 (Particle<number> *p, Particle<number> *q);
	inline number _stacking_energy_n5 (Particle<number> *p, Particle<number> *q);
	inline number _stacking_energy_pq (Particle<number> *p, Particle<number> *q);
	inline number _stacking_energy (Particle<number> *p);
	
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
	CMMC_CPUBackend(IOManager *IO);
	virtual ~CMMC_CPUBackend();

	virtual void get_settings (input_file &inp);
	virtual void print_conf (llint curr_step, bool reduced, bool only_last);
	virtual void print_conf (llint curr_step, bool only_last);
	virtual void print_energy (llint curr_step);
	
	//void init(ifstream &conf_input);
	void init(char conf_filename[256]);

	void sim_step(llint cur_step);
	inline void check_energy();
	inline void check_ops();
	char * get_op_state_str ();
};

#endif /* CMMC_CPUBACKEND_H_ */
