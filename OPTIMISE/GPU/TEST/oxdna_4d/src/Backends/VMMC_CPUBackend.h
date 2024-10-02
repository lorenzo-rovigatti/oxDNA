/**
 * @file    VMMC_CPUBackend.h
 * @date    26/nov/2010
 * @author  flavio
 *
 */

#ifndef VMMC_CPUBACKEND_H_
#define VMMC_CPUBACKEND_H_

#include "MC_CPUBackend.h"
#include "../Utilities/Weights.h"
#include "../Utilities/OrderParameters.h"
#include "../Utilities/Histogram.h"

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)>(b))?(b):(a))

/**
 *
 * @brief This class implements VMMC for the DNA/RNA models
 *
 *
 @verbatim
 [maxclust = <int> (Default: N; maximum number of particles to be moved together. Defaults to the whole system)]
 [small_system = <bool> (Default: false; whether to use an interaction computation suited for small systems.)]
 [preserve_topology = <bool> (Default: false; sets a maximum size for the move attempt to 0.5, which guarantees that the topology of the system is conserved. Also prevents very large moves and might speed up simulations of larger systems, while suppressing diffusion)]
 [umbrella_sampling = <bool> (Default: false; whether to use umbrella sampling)]
 [op_file = <string> (Mandatory if umbrella_sampling is set to true; path to file with the description of the order parameter)]
 [weights_file = <string> (Mandatory if umbrella_sampling is set to true; path to file with the weights to use in umbrella sampling)]
 [last_hist_file = <string> (Optional if umbrella_sampling is set to true, otherwise ignored; Default: last_hist.dat; path to file where the histograms associated with umbrella sampling will be stored. This is printed with the same frequency as the energy file. Should become an observable sooner or later)]
 [traj_hist_file = <string> (Optional if umbrella_sampling is set to true, otherwise ignored; Default: traj_hist.dat; path to file where the series histograms associated with umbrella sampling will be stored, allowing to monitor the time evolution of the histogram and possibly to remove parts of the simulation. This is printed with the same frequency as the energy file. Should become an observable sooner or later)]
 [init_hist_file = <string> (Optional if umbrella_sampling is set to true, otherwise ignored; Default: none; path to a file to load a previous histogram from, useful if one wants to continue a simulation to obtain more statistics.)]
 [extrapolate_hist = <float>,<float>,..., <float> (Optional if umbrella_sampling is set to true, otherwise ignored; Default: none; series of temperatures to which to extrapolate the histograms. They can be given as float in reduced units, or the units can be specified as in the T option)]
 [safe_weights = <bool> (Default: true; whether to check consistency in between order parameter file and weight file. Only used if umbrella_sampling = true)]
 [default_weight = <float> (Default: none; mandatory if safe_weights = true; default weight for states that have no specified weight assigned from the weights file)]
 [skip_hist_zeros = <bool> (Default: false; Wether to skip zero entries in the traj_hist file)]
 [equilibration_steps = <int> (Default: 0; number of steps to ignore to allow for equilibration)]
 @endverbatim
 *
 */
class VMMC_CPUBackend: public MC_CPUBackend {
	struct movestr {
		int seed, type;
		LR_vector t;
		LR_matrix R;
		LR_matrix Rt;
	};

protected:
	bool _have_us;
	bool _reload_hist;
	OrderParameters _op;
	char _op_file[512];
	Weights _w;
	char _weights_file[512];
	Histogram _h;

	bool _safe_weights;
	bool _skip_hist_zeros;
	number _default_weight;

	int _last_move;
	number _U_ext;

	/// potential energy due to stacking
	number _U_stack;

	/// change in stacking potential energy
	number _dU_stack;

	/// change in potential energy
	number _dU;

	inline void store_particle(BaseParticle * src);
	inline void restore_particle(BaseParticle * src);

	int **_neighcells, *_cells, *_vmmc_heads;
	int _vmmc_N_cells, _vmmc_N_cells_side;
	number _vmmc_box_side;

	void _create_cells();
	void _init_cells();
	void _delete_cells();
	inline void _fix_list(int, int, int);
	inline int _get_cell_index(const LR_vector &pos);

	void _print_pos(int);

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

	number ** eijm, **eijm_old;
	bool ** hbijm, **hbijm_old;

	inline number _excluded_volume(const LR_vector &r, number sigma, number rstar, number b, number rc);
	inline number _excluded_volume_faster(const LR_vector &r, const number sigma, const number rstar, const number b, const number rc);

	number _compute_energy_n2();
	void _compute_energy();

	number _particle_particle_bonded_interaction_n5_VMMC(BaseParticle *p, BaseParticle *q, number *stacking_en = 0);

	/**
	 * @brief Computes the bonded interactions given a pair of particles.
	 *
	 * If the third argument (default null) / is set not to be the NULL pointer, the
	 * address it points to is filled with the stacking energy.
	 */
	number _particle_particle_bonded_interaction_n3_VMMC(BaseParticle *p, BaseParticle *q, number *stacking_en = 0);
	number _particle_particle_nonbonded_interaction_VMMC(BaseParticle *p, BaseParticle *q, number *H_energy = 0);

	number build_cluster(movestr *moveptr, int maxsize, int *clust, int *size);
	number build_cluster_cells(movestr *moveptr, int maxsize, int *clust, int *size);
	number build_cluster_small(movestr *moveptr, int maxsize, int *clust, int *size);

	number VMMC_link(double E_new, double E_old) {
		return (1. - exp((1. / this->_T) * (E_old - E_new)));
	}

	inline number _next_rand() {
		return drand48();
	}

	inline void _move_particle(movestr *moveptr, BaseParticle *p, BaseParticle *q);
	void _update_ops();
	void _update_lists();

	inline int cell_neighbours(int myn, int ii) {
		static int x, y, z, nind[3];

		x = myn % _vmmc_N_cells_side;
		y = (myn / _vmmc_N_cells_side) % _vmmc_N_cells_side;
		z = myn / (_vmmc_N_cells_side * _vmmc_N_cells_side);

		nind[0] = (ii / (3 * 3) - 1 + x + _vmmc_N_cells_side) % _vmmc_N_cells_side;
		nind[1] = ((ii / 3) % 3 - 1 + y + _vmmc_N_cells_side) % _vmmc_N_cells_side;
		nind[2] = (ii % 3 - 1 + z + _vmmc_N_cells_side) % _vmmc_N_cells_side;

		return nind[0] + nind[1] * _vmmc_N_cells_side + nind[2] * _vmmc_N_cells_side * _vmmc_N_cells_side;
	}

	std::vector<number> new_en3s, new_stn3s, new_en5s, new_stn5s;

	llint _equilibration_steps;

public:
	VMMC_CPUBackend();
	virtual ~VMMC_CPUBackend();

	virtual void get_settings(input_file &inp);
	virtual void print_conf(bool reduced, bool only_last);
	virtual void print_conf(bool only_last);
	virtual void print_observables();

	void fix_diffusion();

	void init();

	void sim_step();
	inline void check_overlaps();
	inline void check_ops();
	char * get_op_state_str();
};

#endif /* CMMC_CPUBACKEND_H_ */
