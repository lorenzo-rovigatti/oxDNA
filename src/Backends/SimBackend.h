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

#include "../defs.h"
#include "../Observables/ObservableOutput.h"
#include "../Particles/Molecule.h"

#include <cmath>
#include <fstream>
#include <cfloat>
#include <vector>
#include <map>

class IBaseInteraction;
class BaseBox;
class BaseList;
class BaseParticle;
class ObservableOutput;
class ConfigInfo;
class Timer;

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

 [output_prefix = <string> (the name of all output files will be preceded by this prefix, defaults to an empty string)]

 [checkpoint_every = <int> (If > 0, it enables the production of checkpoints, which have a binary format. Beware that trajectories that do have this option enabled will differ from trajectories that do not. If this key is specified, at least one of checkpoint_file and checkpoint_trajectory needs to be specified)]
 [checkpoint_file = <string> (File name for the last checkpoint. If not specified, the last checkpoint will not be printed separately)]
 [checkpoint_trajectory = <string> (File name for the checkpoint trajectory. If not specified, only the last checkpoint will be printed)]
 [reload_from = <string> (checkpoint to reload from. This option is incompatible with the keys conf_file and seed, and requires restart_step_counter=0 as well as binary_initial_conf!=1)]

 @endverbatim
 */
class SimBackend {
protected:
	std::string _backend_info;

	std::shared_ptr<Timer> _mytimer, _obs_timer;

	number _max_io;

	bool _enable_fix_diffusion;

	bool _external_forces;
	std::string _external_filename;

	char _reduced_conf_output_dir[256];

	bool _reseed;

	int _N_strands;
	number _T;
	number _P;
	std::shared_ptr<BaseBox> _box;

	int _conf_interval;
	std::string _conf_filename;
	bool _initial_conf_is_binary;
	bool _back_in_box;
	bool _custom_conf_name;
	char _custom_conf_str[256];
	std::ifstream _conf_input;
	llint _read_conf_step;
	std::string _checkpoint_file;
	std::string _checkpoint_traj;
	bool _restart_step_counter;

	/// Vector of ObservableOutput used to manage the simulation output
	std::vector<ObservableOutputPtr> _obs_outputs;
	ObservableOutputPtr _obs_output_stdout;
	ObservableOutputPtr _obs_output_file;
	ObservableOutputPtr _obs_output_trajectory;
	ObservableOutputPtr _obs_output_last_conf;
	ObservableOutputPtr _obs_output_last_conf_bin;
	ObservableOutputPtr _obs_output_reduced_conf;
	ObservableOutputPtr _obs_output_checkpoints;
	ObservableOutputPtr _obs_output_last_checkpoint;

	/// Shared pointer to the interaction manager
	InteractionPtr _interaction;

	/// Pointer to the list manager
	ListPtr _lists;

	number _rcut;
	number _sqr_rcut;

	/// potential energy
	number _U;

	/// array of pointers to particle objects
	std::vector<BaseParticle *> _particles;

	std::vector<std::shared_ptr<Molecule>> _molecules;

	/// object that stores pointers to a few important variables that need to be shared with other objects
	ConfigInfo *_config_info = nullptr;

	int _N_updates;
	int _confs_to_skip;
	llint _bytes_to_skip;

	/**
	 * Reads from _conf_input three numbers in ascii or binary format, depending on the
	 * value of the parameter.
	 * @param binary whether _conf_input has been open in ascii or binary format
	 * @return a vector containing the three numbers read
	 */
	LR_vector _read_next_binary_vector();

	virtual void _on_T_update();

public:
	SimBackend();
	virtual ~SimBackend();

	/**
	 * @brief Loads all the required settings from the input file.
	 *
	 * @param inp
	 */
	virtual void get_settings(input_file &inp);
	/**
	 * @brief Initializes the backend.
	 *
	 * @param path
	 */
	virtual void init();

	/**
	 * @brief Reads the next configuration from the conf_file.
	 *
	 * It can read it either in binary or ascii format.
	 *
	 * @param binary whether conf_file is to be parsed in ascii or binary format
	 * @return true if the operation was successful, false otherwise
	 */
	virtual bool read_next_configuration(bool binary=false);

	int N() {
		return _particles.size();
	}

	/**
	 * @brief Returns how many times verlet lists have been updated.
	 *
	 * @return number of Verlet lists updates
	 */
	int get_N_updates() {
		return _N_updates;
	}

	virtual void fix_diffusion();
	virtual void print_equilibration_info();

	void add_output(ObservableOutputPtr new_output);
	void remove_output(std::string output_file);

	/**
	 * @brief Prints the observables attached to the backend.
	 */
	virtual void print_observables();

	virtual void update_observables_data();

	virtual void print_conf(bool reduced=false, bool only_last=false);

	/**
	 * @brief Performs a simulation step.
	 */
	virtual void sim_step() = 0;

	/**
	 * @brief Synchronize the simulation data with the data structures that are used to analyse/print the current simulation status.
	 */
	virtual void apply_simulation_data_changes();

	/**
	 * @brief Update the simulation data, so that changes done to the data structures are taken into account by the simulation engine.
	 */
	virtual void apply_changes_to_simulation_data();

	long long int current_step() {
		return _config_info->curr_step;
	}

	void increment_current_step() {
		_config_info->curr_step++;
	}

	llint start_step_from_file;
};

#endif /* SIMBACKEND_H_ */
