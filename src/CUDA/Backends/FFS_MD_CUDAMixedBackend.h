/*
 * FFS_MD_CUDAMixedBackend.h
 *
 *  Created on: 18/apr/2013
 *      Author: Ben Snodin
 */

#ifndef FFS_MD_CUDAMIXEDBACKEND_H_
#define FFS_MD_CUDAMIXEDBACKEND_H_

#define FORWARD 1
#define BACKWARD 0

#include "MD_CUDAMixedBackend.h"
#include "../../Utilities/OrderParameters.h"
#include "../../Backends/FFS_MD_CPUBackend.h" // for parsed_condition and parsed_expression defs

/* This struct holds all the data (mostly host and device pointers) for a set of
 * conditions. We use a struct to do this because there might be different sets
 * of conditions e.g. for flux generation there is a set of forward conditions 
 * and a set of backward conditions; additionally each master condition element 
 * will have its own SimpleConditions struct containing to the conditions it is 
 * responsible for.
 */
struct SimpleConditions {
	// string containing name of condition set
	char type[256];

	// variables needed to record host 2d variable-width arrays for hb_cond, 
	// nearhb_cond and dist_cond
	int hb_cond_len;
	int *h_hb_cond_lens;
	int *h_hb_cond_rows;
	int nearhb_cond_len;
	int *h_nearhb_cond_lens;
	int *h_nearhb_cond_rows;
	int dist_cond_len;
	int *h_dist_cond_lens;
	int *h_dist_cond_rows;

	// host arrays for magnitudes and types of condition for hb_cond, 
	// nearhb_cond and dist_cond
	float *h_dist_cond_mags;
	int *h_dist_cond_types;
	float *h_hb_cond_mags;
	int *h_hb_cond_types;
	float *h_nearhb_cond_mags;
	int *h_nearhb_cond_types;

	// variables needed to record device 2d variable-width arrays for hb_cond, 
	// nearhb_cond and dist_cond
	int *d_dist_cond_lens;
	int *d_dist_cond_rows;
	int *d_hb_cond_lens;
	int *d_hb_cond_rows;
	int *d_nearhb_cond_lens;
	int *d_nearhb_cond_rows;

	// device arrays for magnitudes and types of condition for hb_cond, 
	// nearhb_cond and dist_cond
	float *d_hb_cond_mags;
	int *d_hb_cond_types;
	float *d_nearhb_cond_mags;
	int *d_nearhb_cond_types;
	float *d_dist_cond_mags;
	int *d_dist_cond_types;

	// each set holds the indices of all the conditions of that type
	// stop_set_count is the total c_number of stopping conditions
	int *nearhb_cond_set;
	int *hb_cond_set;
	int *dist_cond_set;
	int stop_set_count;

	// pointers to send stopping condition status from device to host
	// stop_length gives the length of the corresponding array
	bool *h_ffs_stop;
	bool *d_ffs_stop;
	int stop_length;
};

/*
 * Each line in the conditions file within a master condition definition is
 * stored in a master_condition_element. NB a single master condition may
 * contain elements that check for different condition types (i.e. their names
 * have different prefixes). e.g.
 *
 * master1 = {
 * condition_set_a >= 3
 * condition_set_b < 10
 * }
 *
 */
struct master_condition_element {
	// simple conditions that this master condition will check
	SimpleConditions simple_conditions;

	// conditions_name, value and type define the master condition e.g. 
	// "melting > 10" corresponds to prefix=melting type=">" value = 10
	// and means more than 10 conditions with name prefix 'melting' must be
	// satisfied
	char conditions_name[256];
	int value;
	int type;

	bool parse_master_element(const char *expression, const char *fname);
};

/* A master condition is a set of conditions based on the collective status of other
 * conditions. For example it is useful if you have 100 possible conditions that
 * might be satisfied, and you want to stop the simulation if any 10 of those
 * conditions are satisfied. e.g.
 *
 * master1 = {
 * my_condition >= 10
 * }
 *
 * my_condition1 = {
 * dist1 > 5
 * }
 * my_condition2 = {
 * dist2 > 5
 * }
 * my_condition3 = {
 * dist3 > 5
 * }
 * ...
 */
struct master_condition {
	std::vector<master_condition_element> elements;
	char name[256];
	bool parse_master_condition(const char *expression, const char *fname);
};

/**
 * @brief CUDA backend with forward flux sampling capability
 *
 * This class is derived from the CUDA MD mixed backend with 
 * added forward flux sampling capability through a modified
 * sim_step() method and a few others. It can run both flux 
 * calculation and shooting runs. It works with mixed precision
 * on CUDA only.
 *
 * Input options:
 * 
 * @verbatim
 backend = CUDA (For CUDA FFS -- NB unlike the CPU implementation, the CUDA implementation does not print extra columns with the current order parameter values whenever the energy is printed)
 sim_type = FFS_MD (This must be set for an FFS simulation)
 order_parameters_file = <string> (path to the order parameters file)
 ffs_file = <string> (path to the file with the simulation stopping conditions. Optionally, one may use 'master conditions' (CUDA FFS only), which allow one to more easily handle very high dimensional order parameters. See the EXAMPLES/CUDA_FFS/README file for more information)
 [ffs_generate_flux = <bool> (CUDA FFS only. Default: False; if False, the simulation will run until a stopping condition is reached; if True, a flux generation simulation will be run, in which case reaching a condition will cause a configuration to be saved but will not terminate the simulation. In the stopping condition file, the conditions must be labelled forward1, forward2, ... (for the forward conditions); and backward1, backward2, ... (for the backward conditions), ... instead of condition1, condition2, ... . To get standard flux generation, set the forward and backward conditions to correspond to crossing the same interface (and use conditions corresponding to different interfaces for Tom's flux generation). As with the single shooting run mode, the name of the condition crossed will be printed to stderr each time.)]
 [gen_flux_save_every = <integer> (CUDA FFS only. Mandatory if ffs_generate_flux is True; save a configuration for 1 in every N forward crossings)]
 [gen_flux_total_crossings = <integer> (CUDA FFS only. Mandatory if ffs_generate_flux is True; stop the simulation after N crossings achieved)]
 [gen_flux_conf_prefix = <string> (CUDA FFS only. Mandatory if ffs_generate_flux is True; the prefix used for the file names of configurations corresponding to the saved forward crossings. Counting starts at zero so the 3rd crossing configuration will be saved as MY_PREFIX_N2.dat)]
 [gen_flux_debug = <bool> (CUDA FFS only. Default: False; In a flux generation simulation, set to true to save backward-crossing configurations for debugging)]
 [check_initial_state = <bool> (CUDA FFS only. Default: False; in a flux generation simulation, set to true to turn on initial state checking. In this mode an initial configuration that crosses the forward conditions after only 1 step will cause the code to complain and exit. Useful for checking that a flux generation simulation does not start out of the A-state)]
 [die_on_unexpected_master = <bool> (CUDA FFS only. Default: False; in a flux generation simulation that uses master conditions, set to true to cause the simulation to die if any master conditions except master_forward1 or master_backward1 are reached. Useful for checking that a flux generation simulation does not enter any unwanted free energy basins (i.e. other than the initial state and the desired final state))]
 [unexpected_master_prefix = <string> (CUDA FFS only. Mandatory if die_on_unexpected_master is True; the prefix used for the file names of configurations corresponding to reaching any unexpected master conditions (see die_on_unexpected_master).)]
 @endverbatim
 */
class FFS_MD_CUDAMixedBackend: public CUDAMixedBackend {
protected:
	// these are parameters read from the input file
	char _order_parameters_file[256];
	char _ffs_file[256];
	bool _gen_flux;
	int _gen_flux_save_every;
	int _gen_flux_desired_cc;
	int _gen_flux_debug;
	char _conf_prefix[256];
	int _check_initial_state;
	int _die_on_unexpected_master;
	bool _unexpected_master_reached;
	char _unexpected_master_name[256];
	char _unexpected_master_prefix[256];
	char _state_str[2048];
	int _print_energy_every;

	int _gen_flux_cross_count;
	int _gen_flux_saved_cross_count;
	int _flux_direction;
	bool _use_master_conditions;
	OrderParameters _op;

	// conditions
	SimpleConditions _sc;
	SimpleConditions _sc_fwd;
	SimpleConditions _sc_bwd;

	ObservableOutputPtr _obs_output_custom_conf;

	// for init_from_ffs_file
	std::vector<parsed_condition> _conditions;
	std::vector<parsed_condition> _forward_conditions;
	std::vector<parsed_condition> _backward_conditions;
	std::vector<master_condition> _master_conditions;
	std::vector<master_condition> _master_forward_conditions;
	std::vector<master_condition> _master_backward_conditions;

	CUDA_kernel_cfg _ffs_hb_precalc_kernel_cfg;
	CUDA_kernel_cfg _ffs_dist_precalc_kernel_cfg;
	CUDA_kernel_cfg _ffs_hb_eval_kernel_cfg;
	CUDA_kernel_cfg _ffs_dist_eval_kernel_cfg;

	// containers for the order parameters
	int _n_hb_pairs;
	int _n_hb_regions;
	int _n_dist_pairs;
	int _n_dist_regions;
	float *_h_hb_energies;
	bool *_h_nearhb_states;
	float *_h_op_dists;
	int *_d_hb_pairs1;
	int *_d_hb_pairs2;
	int *_d_dist_pairs1;
	int *_d_dist_pairs2;
	float *_d_op_dists;
	float *_d_hb_energies;
	float *_d_hb_cutoffs;
	bool *_h_region_is_nearhb;

	int *_d_dist_region_lens;
	int *_d_dist_region_rows;
	int *_d_hb_region_lens;
	int *_d_hb_region_rows;
	bool *_d_region_is_nearhb;
	bool *_d_nearhb_states;

	/**
	 * @brief Parse the conditions file
	 *
	 * @param fname name of the conditions file
	 */
	void _init_ffs_from_file(const char *fname);

	/**
	 * @brief Configure the kernel used to check stopping conditions
	 *
	 * @param kernel_cfg pointer to CUDA kernel config
	 * @param total_threads total c_number of threads required
	 */
	void _init_ffs_kernel_config(CUDA_kernel_cfg *kernel_cfg, int total_threads);

	/**
	 * @brief Make sure the required CUDA constants are available to the kernels that run the simulation and check the stopping conditions
	 */
	void _init_CUDA_MD_symbols();

	/**
	 * @brief Return a pointer to an array of indices for the beginning of each row for a variable-width 2D array A using an array storing the lengths of each row of A
	 *
	 * @param rows_len c_number of rows
	 * @param lens array containing length of each row
	 */
	int *_get_2D_rows(int rows_len, int *lens);

	/**
	 * @brief Parse conditions with a certain prefix -- overloaded for a vector of parsed_condition pointers (parse a set of simple conditions)
	 *
	 * @param fname name of the conditions file
	 * @param condition_set_type name of the condition prefix to be parsed
	 * @param conditions a vector to store the conditions
	 */
	bool _read_conditions(const char *fname, const char *condition_set_type, std::vector<parsed_condition> *conditions);

	/**
	 * @brief Parse conditions with a certain prefix -- overloaded for a vector of master_condition pointers (parse a set of master conditions)
	 *
	 * @param fname name of the conditions file
	 * @param condition_set_type name of the condition prefix to be parsed
	 * @param conditions a vector to store the conditions
	 */
	bool _read_conditions(const char *fname, const char *condition_set_type, std::vector<master_condition> *conditions);

	/**
	 * @brief Copy the data for a particular set of conditions to the device and prepare the host structs. The nice
	 * C++ vector of structs that is built by read_conditions() is converted to a struct of arrays, which are then copied
	 * to the device.
	 *
	 * @param conditions set of conditions created by read_conditions()
	 * @param type the prefix for the name of this set of conditions
	 */
	SimpleConditions _get_simple_conditions(std::vector<parsed_condition> conditions, const char type[256]);

	/**
	 * @brief For each master condition, prepare the SimpleConditions struct for each master condition element
	 *
	 * @param master_conditions all the master conditions
	 * @param fname name of the conditions file
	 */
	void _master_conditions_prepare_simple_conditions(vector<master_condition> *master_conditions, const char *fname);

	/**
	 * @brief For a given set of simple conditions, return the total c_number of conditions satisfied. Useful for master condition element evaluation
	 *
	 * @param sc contains the simple conditions to check
	 * @param suppress_logging optionally don't print to stdout when condition is satisfied (currently on for flux generation, off for master condition evaluation)
	 */
	int _test_crossing(SimpleConditions sc, bool suppress_logging = false);

	/**
	 * @brief Launch kernels to check the current order parameter states
	 */
	void _eval_order_parameter_states();

	/**
	 * @brief Check the status of the stop conditions for a set of simple conditions on the device and copy the result to a host array
	 *
	 * @param sc set of conditions
	 */
	void _eval_stop_conditions(SimpleConditions sc);

	/**
	 * @brief Get ready to save a configuration by copying data from the device and setting the configuration name
	 *
	 * @param conf_name configuration name
	 */
	void _prepare_configuration(char *conf_name);

	/**
	 * @brief Test all master conditions and return True and print to the screen if one of them is met. Also check for
	 * unexpected master crossings and set a flag if found.
	 *
	 * @param master_conditions master conditions to test
	 */
	bool _test_master_conditions(vector<master_condition> master_conditions);

	/**
	 * @brief For each element of a particular master condition, check the state of all of the element's condition sets.
	 * If the right c_number of each element's conditions are satisfied, return true, otherwise return false
	 *
	 * @param master_condition the master condition to test
	 */
	bool _test_master_condition(master_condition master_condition);

	/**
	 * @brief Print to the screen the names of all simple conditions corresponding to a certain master condition that has been reached
	 *
	 * @param master_condition the master condition
	 */
	void _log_master_state(master_condition master_condition);

	/**
	 * @brief Print information to screen and save a configuration. To be used when an unexpected master condition is reached.
	 */
	void _handle_unexpected_master();

	/**
	 * @brief Free all memory allocation for master_conditions
	 *
	 * @brief master_conditions master conditions to free
	 */
	void _free_master_conditions(vector<master_condition> master_conditions);

	/**
	 * @brief Free all memory allocation for a SimpleConditions struct
	 *
	 * @param sc struct to be freed
	 */
	void _free_simple_conditions(SimpleConditions sc);

	/**
	 * @brief Return true if the stopping conditions have been met and the simulation should be 
	 * terminated. Run every step, conditions are computed on the device
	 */
	bool _check_stop();

public:
	/**
	 * @brief Constructor. Sets a load of default values
	 */
	FFS_MD_CUDAMixedBackend();
	virtual ~FFS_MD_CUDAMixedBackend();

	/**
	 * @brief Get the order parameters from a file and copy them to the device.
	 *
	 * @param conf_filename order parameters file
	 */
	void init();

	/**
	 * @brief Read settings from the input file
	 *
	 * @param inp input file
	 */
	void get_settings(input_file &inp);

	/**
	 * @brief Move the simulation forward by one step; check for a stopping condition with _check_stop() 
	 * and end the simulation if it returns True.
	 * 
	 * @param curr_step current step
	 */
	void sim_step();

	/**
	 * @brief Get a string showing the current state of the order parameters
	 */
	char * get_op_state_str(void);

	/**
	 * @brief Prints observables; makes sure order parameter values are up-to-date
	 *
	 * @param curr_step current step
	 */
	virtual void print_observables();

	/**
	 * @brief Print a list of order parameter names and values to the string str. Designed to reproduce
	 * the output at the end of a CPU FFS simulation 31/12/14
	 *
	 * @param str string to be written to, memory should already be allocated
	 */
	void sprintf_names_and_values(char *str);
};

#endif /* FFS_MD_CUDAMIXEDBACKEND_H_ */
