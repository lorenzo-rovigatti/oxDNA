/**
 * @file    SimBackend.h
 * @date    03/set/2010
 * @author  lorenzo
 *
 *
 */

#ifndef SIMBACKEND_H_
#define SIMBACKEND_H_

#define SIM_MD 0
#define SIM_MC 1

#include <cmath>
#include <fstream>
#include <cfloat>
#include <vector>

#include "../defs.h"
#include "../Particles/BaseParticle.h"
#include "../Observables/ObservableOutput.h"
#include "../Lists/BaseList.h"
#include "../Boxes/BaseBox.h"
#include "../Utilities/Timings.h"

using namespace std;

template <typename number> class IBaseInteraction;

/**
 * @brief Defines the backend interface. It is an abstract class.
 */
class ISimBackend {
public:
	ISimBackend() : _start_step_from_file(0) {};
	virtual ~ISimBackend() {};

	llint _start_step_from_file;

	/**
	 * @brief Loads all the required settings from the input file.
	 *
	 * @param inp
	 */
	virtual void get_settings(input_file &inp) = 0;

	/**
	 * @brief Initializes the backend.
	 *
	 * @param path
	 */
	virtual void init() = 0;

	/**
	 * @brief Performs a simulation step.
	 *
	 * @param curr_step
	 */
	virtual void sim_step(llint curr_step) = 0;
	virtual void print_conf(llint curr_step, bool reduced=false, bool only_last=false) = 0;

	/**
	 * @brief Returns how many times verlet lists have been updated.
	 *
	 * @return number of Verlet lists updates
	 */
	virtual int get_N_updates() = 0;

	/**
	 * @brief Prints the observables attached to the backend.
	 *
	 * @param curr_step
	 */
	virtual void print_observables(llint curr_step) = 0;

	virtual void fix_diffusion() = 0;
};

/**
 * @brief Base backend class. Every backend inherits from this class.
 *
 * @verbatim
T = <float> (temperature of the simulation. It can be expressed in simulation units or kelvin (append a k or K after the value) or celsius (append a c or C after the value).)

[fix_diffusion = <bool> (if true, particles that leave the simulation box are brought back in via periodic boundary conditions. Defaults to true.)]
[seed = <int> (seed for the random number generator. On Unix systems, defaults to either a number from /dev/urandom or to time(NULL))]

[confs_to_skip = <int> (how many configurations should be skipped before using the next one as the initial configuration, defaults to 0)]
restart_step_counter = <boolean> (false means that the step counter will start from the value read in the configuration file, true means that the step counter will start from 0)

[external_forces = <bool> (specifies whether there are external forces acting on the nucleotides or not. If it is set to 1, then a file which specifies the external forces' configuration has to be provided (see external_forces_file))]
[external_forces_file = <path> (specifies the file containing all the external forces' configurations. Currently there are six supported force types: string, twist, trap, repulsion_plane, repulsion_plane_moving and mutual_trap (see EXAMPLES/TRAPS for some examples))]

[back_in_box = <bool> (whether particles should be brought back into the box when a configuration is printed or not, defaults to false)]

[lastconf_file = <path> (path to the file where the last configuration will be dumped)]
trajectory_file = <path> (path to the file which will contain the output trajectory of the simulation)

[binary_initial_conf = <bool> (whether the initial configuration is a binary configuration or not, defaults to false)]
[lastconf_file_bin = <path> (path to the file where the last configuration will be printed in binary format, if not specified no binary configurations will be printed)]

[print_reduced_conf_every = <int> (every how many time steps configurations containing only the centres of mass of the strands should be printed. If 0, no reduced configurations will be printed)]
[reduced_conf_output_dir = <path> (path to the folder where reduced configurations will be printed)]

[no_stdout_energy = <bool> (if true oxDNA will not print the default simulation output, including the energy, to stdout. Defaults to false)]

[print_timings = <bool> (whether oxDNA should print out to a file performance timings at the end of the simulation or not, defaults to false)]
[timings_filename = <path> (path to the file where timings will be printed)]

[output_prefix = <string> (the name of all output files will be preceded by this prefix, defaults to an empty string)]

[checkpoint_every = <int> (If > 0, it enables the production of checkpoints, which have a binary format. Beware that trajectories that do have this option enabled will differ from trajectories that do not. If this key is specified, at least one of checkpoint_file and checkpoint_trajectory needs to be specified)]
[checkpoint_file = <string> (File name for the last checkpoint. If not specified, the last checkpoint will not be printed separately)]
[checkpoint_trajectory = <string> (File name for the checkpoint trajectory. If not specified, only the last checkpoint will be printed)]
[reload_from = <string> (checkpoint to reload from. This option is incompatible with the keys conf_file and seed, and requires restart_step_counter=0 as well as binary_initial_conf!=1)]

@endverbatim
 */
template<typename number>
class SimBackend: public ISimBackend{
protected:
	std::string _backend_info;

	Timer * _mytimer;

	number _max_io;

	bool _enable_fix_diffusion;
	bool _print_timings;
	char _timings_filename[256];
	bool _is_CUDA_sim;

	bool _external_forces;
	char _external_filename[256];

	char _reduced_conf_output_dir[256];

	int _sim_type;

	bool _reseed;

	int _N;
	int _N_strands;
	number _T;
	number _P;
	number _box_side;
	BaseBox<number> *_box;

	int _conf_interval;
	std::string _conf_filename;
	bool _initial_conf_is_binary;
	bool _back_in_box;
	bool _custom_conf_name;
	char _custom_conf_str[256];
	ifstream _conf_input;
	llint _read_conf_step;
	std::string _checkpoint_file;
	std::string _checkpoint_traj;

	/// Vector of ObservableOutput used to manage the simulation output
	vector<ObservableOutput<number> *> _obs_outputs;
	ObservableOutput<number> *_obs_output_stdout;
	ObservableOutput<number> *_obs_output_file;
	ObservableOutput<number> *_obs_output_trajectory;
	ObservableOutput<number> *_obs_output_last_conf;
	ObservableOutput<number> *_obs_output_last_conf_bin;
	ObservableOutput<number> *_obs_output_reduced_conf;
	ObservableOutput<number> *_obs_output_checkpoints;
	ObservableOutput<number> *_obs_output_last_checkpoint;

	/// Pointer to the interaction manager
	IBaseInteraction<number> *_interaction;

	/// Pointer to the list manager
	BaseList<number> *_lists;

	number _rcut;
	number _sqr_rcut;

	/// potential energy
	number _U;

	/// kinetic energy
	number _K;

	/// potential energy due to hydrogen bonding
	number _U_hydr;

	/// potential energy due to stacking
	number _U_stack;

	/// change in stacking potential energy
	number _dU_stack;

	/// change in potential energy
	number _dU;

	/// array of pointers to particle objects
	BaseParticle<number> **_particles;

	/// object that stores pointers to a few important variables that need to be shared with other objects
	ConfigInfo<number> _config_info;

	void _get_number_settings(input_file &inp);

	int _N_updates;
	int _confs_to_skip;

	void _read_external_forces();

	/**
	 * Reads from _conf_input three numbers in ascii or binary format, depending on the
	 * value of the parameter.
	 * @param binary whether _conf_input has been open in ascii or binary format
	 * @return a vector containing the three numbers read
	 */
	template <typename n_number> LR_vector<n_number> _read_next_vector(bool binary);

	/**
	 * @brief Reads the next configuration from the conf_file.
	 *
	 * It can read it either in binary or ascii format.
	 *
	 * @param binary whether conf_file is to be parsed in ascii or binary format
	 * @return true if the operation was successful, false otherwise
	 */
	bool _read_next_configuration(bool binary=false);

	int _get_N_from_conf(ifstream &conf_input);

	/**
	 * @brief Prints all the observables that are ready, i.e. whose is_ready(curr_step) method returns true.
	 *
	 * @param curr_step
	 */
	virtual void _print_ready_observables(llint curr_step);

public:
	SimBackend();
	virtual ~SimBackend();

	virtual void get_settings(input_file &inp);
	virtual void init();

	int get_N_updates() {return _N_updates; }
	virtual void fix_diffusion();
	virtual void print_observables(llint curr_step);
	virtual void print_conf(llint curr_step, bool reduced=false, bool only_last=false);
};

#endif /* SIMBACKEND_H_ */
