/*
 * CMMC_CPUBackend.h
 *
 *  Created on: 26/nov/2010
 *      Author: flavio 
 */

#ifndef VMMC_CPUBACKEND_H_
#define VMMC_CPUBACKEND_H_

#include "MCBackend.h"
#include "MC_CPUBackend.h"
#include "Weights.h"
#include "OrderParameters.h"
#include "Histogram.h"

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)>(b))?(b):(a))

template<typename number>
struct movestr {
    int seed, type;
    LR_vector<number> t;
    LR_matrix<number> R;
    LR_matrix<number> Rt;
};

template<typename number>
class VMMC_CPUBackend: public MC_CPUBackend<number> {
protected:
	Particle<number> *_particles_old;

	bool _have_us;
	bool _reload_hist;
	OrderParameters _op;
	char _op_file[512];
	Weights _w;
	char _weights_file[512];
	Histogram _h;

	int _last_move;

	//void store_particle (Particle<number> *, Particle<number> *);
	inline void store_particle (Particle<number> * src);
	inline void restore_particle (Particle<number> * src);
	
	/*
	int * _heads, ** _neighcells, * _cells, _N_cells, _N_cells_side;
	*/
	int **_neighcells, * _cells, * _vmmc_heads;
	int _vmmc_N_cells, _vmmc_N_cells_side;
	void _create_cells();
	void _init_cells();
	void _delete_cells();
	inline void _fix_list(int, int, int);
	inline int _get_cell_index(const LR_vector<number> &pos);

	void _print_pos (int);

	char _traj_hist_file[512];
	char _last_hist_file[512];
	char _init_hist_file[512];
	char _state_str[512];

	int _maxclust;

	bool _preserve_topology, _small_system;
	number _max_move_size, _max_move_size_sqr;
	
	bool _just_updated_lists;
	bool _reject_prelinks;
	int _netemps;
	double * _etemps;

	number ** eijm, ** eijm_old;
	bool ** hbijm, ** hbijm_old;
	
	inline number _excluded_volume(const LR_vector<number> &r, number sigma, number rstar, number b, number rc);
	inline number _excluded_volume_faster(const LR_vector<number> &r, const number sigma, const number rstar, const number b, const number rc);

	number _compute_energy_n2();
	void _compute_energy();

	inline number _particle_particle_bonded_interaction_n5 (Particle<number> *p, Particle<number> *q,number *stacking_en = 0);
	inline number _particle_particle_bonded_interaction_n3 (Particle<number> *p, Particle<number> *q,number *stacking_en = 0);
	inline number _particle_particle_interaction (Particle<number> *p, Particle<number> *q, number *H_energy=0);
	
	number build_cluster (movestr<number> * moveptr, int maxsize, int * clust, int * size);
	number build_cluster_cells (movestr<number> * moveptr, int maxsize, int * clust, int * size);
	number build_cluster_small (movestr<number> * moveptr, int maxsize, int * clust, int * size);
	number build_cluster_smallish (movestr<number> * moveptr, int maxsize, int * clust, int * size);

	number VMMC_link(double E_new, double E_old) { return (1. - exp((1. / this->_T) * (E_old - E_new)));}

	inline number _next_rand (void) const { return drand48(); }

	inline void _move_particle(movestr<number> * moveptr, Particle<number> *p);
	//void _r_move_particle(movestr<number> * moveptr, Particle<number> *p);
	inline void _fill_e_neigh (Particle<number> *, Particle<number> *, number, int);
	inline void _fill_h_bonds (Particle<number> *, Particle<number> *, bool);
	void _update_metainfo ();
	void _check_metainfo ();
	void _check_old_metainfo ();
	void _update_ops ();
	void _update_lists ();

	number * new_en3s, * new_stn3s, * new_en5s, * new_stn5s;

public:
	VMMC_CPUBackend(IOManager *IO);
	virtual ~VMMC_CPUBackend();

	virtual void get_settings (input_file &inp);
	virtual void print_conf (llint curr_step, bool reduced, bool only_last);
	virtual void print_conf (llint curr_step, bool only_last);
	virtual void print_energy (llint curr_step);
	
	//void init(ifstream &conf_input);
	void init(char conf_filename[256]);

	void sim_step(llint cur_step);
	inline void check_energy();
	inline void check_overlaps();
	inline void check_ops();
	char * get_op_state_str ();
};

#endif /* CMMC_CPUBACKEND_H_ */
