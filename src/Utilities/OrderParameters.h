/**
 * @file    OrderParameters.h
 * @date    Feb 10, 2012
 * @author  sulc
 */

#ifndef ORDERPARAMETERS_H_
#define ORDERPARAMETERS_H_

#include "../Utilities/Utils.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"
#include "../Particles/DNANucleotide.h"
#include "../Boxes/BaseBox.h"

#include <set>
#include <utility>
#include <vector>
#include <string>
#include <map>
#include <list>
#include <iterator>

using namespace std;

typedef pair<int, int> base_pair;
typedef vector<base_pair> vector_of_pairs;

struct classcomp {
	bool operator()(const base_pair& a, const base_pair& b) const {
		if(a.first < b.first)
		return true;
		else if(a.first == b.first && a.second < b.second)
		return true;
		else
		return false;
	}
};

typedef std::set<base_pair, classcomp> set_of_pairs;

struct HBParameter {
	/// this is what defines order parameter: it is loaded from external init file
	set_of_pairs counted_pairs;
	int current_value;
	int stored_value;
	std::string name;

	double cutoff;

	HBParameter();

	void save_pair_as_parameter(int a, int b);

	void set_current_value(int new_value) {
		current_value = new_value;
	}
	int store_current_value() {
		stored_value = current_value;
		return current_value;
	}
	int restore_current_value() {
		current_value = stored_value;
		return current_value;
	}

	double get_cutoff() {
		return cutoff;
	}
	void set_cutoff(double newc) {
		cutoff = newc;
	}
	int get_state() {
		return current_value;
	}

	int get_max_state() {
		return counted_pairs.size();
	}

	std::string get_name() {
		return name;
	}
	void set_name(std::string newname) {
		name = newname;
	}

	int plus_pair(base_pair& bp, double energy = -16);
	int minus_pair(base_pair& bp);
	void reset() {
		current_value = 0;
	}

};

/**
 * @brief op that saves minimum distance of base distances between pairs
 */
struct MinDistanceParameter {
	/// this is what defines order parameter: it is loaded from external init file
	vector_of_pairs counted_pairs;
	double current_value;
	double stored_value;
	std::string name;

	// the sub-type is used to decide what to actually measure 
	int _sub_type;

	//flag for use of COM-COM distances instead of base-base distances:
	bool _use_COM;

	// added to discretize it
	/// index of the state
	int state_index;

	/// to be used in the store/restore functions
	int stored_state_index;

	/// number of different states
	int n_states;

	/// interfaces between states
	vector<double> interfaces;

	MinDistanceParameter();
	void save_pair_as_parameter(int a, int b);

	void reset() {
		current_value = 1.;
		state_index = -1;
	}

	void set_current_value(int new_value) {
		current_value = new_value;
	}

	double store_current_value() {
		stored_value = current_value;
		stored_state_index = state_index;
		return current_value;
	}

	double restore_current_value() {
		state_index = stored_state_index;
		current_value = stored_value;
		return current_value;
	}

	double get_state() {
		return current_value;
	}

	int get_state_index() {
		return state_index;
	}

	std::string get_name() {
		return name;
	}
	void set_name(std::string newname) {
		name = newname;
	}

	int calculate_state(std::vector<BaseParticle *> &particle_list, BaseBox * box);
};

/**
 * Example of order_parameter file:
 * {
 * order_parameter = bond
 * name = all_bonds_01
 * pair1 = 14, 20
 * pair2 = 14, 16
 * pair3 = 17, 41
 * }
 * {
 * order_parameter = mindistance
 * name = dist_strand_12
 * pair1 = 14,58
 * pair2 = 12,45
 * }
 *
 */

class OrderParameters {
protected:
	std::map<std::string, int> _hb_parnames;
	std::map<std::string, int> _dist_parnames;

	vector<HBParameter> _hb_parameters;
	int _hb_parameters_count;
	int *_hb_states;

	vector<MinDistanceParameter> _distance_parameters;
	int _distance_parameters_count;
	double *_distance_states;

	int _all_states_count;
	int * _all_states;

	int _log_level;

public:
	OrderParameters();

	///access functions
	int get_hb_parameters_count() {
		return _hb_parameters_count;
	}

	int get_distance_parameters_count() {
		return _distance_parameters_count;
	}

	int get_all_parameters_count() {
		//return _distance_parameters_count + _hb_parameters_count;
		return _all_states_count;
	}

	// functions to extract a list of particle pair ids that are involved in either a hydrogen
	// bond or a minimum distance order parameter
	void get_dist_pairs(int *out1, int *out2) {
		int ii = 0;
		int op_count = get_distance_parameters_count();
		for(int op_ind = 0; op_ind < op_count; op_ind++) {
			for(vector_of_pairs::iterator i = _distance_parameters[op_ind].counted_pairs.begin();
					i != _distance_parameters[op_ind].counted_pairs.end(); i++) {
				out1[ii] = (*i).first;
				out2[ii] = (*i).second;
				ii++;
			}
		}
	}

	void get_hb_pairs(int *out1, int *out2) {
		int ii = 0;
		int op_count = get_hb_parameters_count();
		for(int op_ind = 0; op_ind < op_count; op_ind++) {
			for(set_of_pairs::iterator i = _hb_parameters[op_ind].counted_pairs.begin();
					i != _hb_parameters[op_ind].counted_pairs.end(); i++) {
				out1[ii] = (*i).first;
				out2[ii] = (*i).second;
				ii++;
			}
		}
	}

	void get_dist_pairs_count(int *counts) {
		int op_count = get_distance_parameters_count();
		for(int op_ind = 0; op_ind < op_count; op_ind++) {
			counts[op_ind] = 0;
			for(vector_of_pairs::iterator i = _distance_parameters[op_ind].counted_pairs.begin();
					i != _distance_parameters[op_ind].counted_pairs.end(); i++) {
				counts[op_ind]++;
			}
		}
	}

	void get_hb_pairs_count(int *counts) {
		int op_count = get_hb_parameters_count();
		for(int op_ind = 0; op_ind < op_count; op_ind++) {
			counts[op_ind] = 0;
			for(set_of_pairs::iterator i = _hb_parameters[op_ind].counted_pairs.begin();
					i != _hb_parameters[op_ind].counted_pairs.end(); i++) {
				counts[op_ind]++;
			}
		}
	}

	int get_dist_pairs_total_count() {
		int ii = 0;
		int op_count = get_distance_parameters_count();
		for(int op_ind = 0; op_ind < op_count; op_ind++) {
			for(vector_of_pairs::iterator i = _distance_parameters[op_ind].counted_pairs.begin();
					i != _distance_parameters[op_ind].counted_pairs.end(); i++) {
				ii++;
			}
		}
		return ii;
	}

	int get_hb_pairs_total_count() {
		int ii = 0;
		int op_count = get_hb_parameters_count();
		for(int op_ind = 0; op_ind < op_count; op_ind++) {
			for(set_of_pairs::iterator i = _hb_parameters[op_ind].counted_pairs.begin();
					i != _hb_parameters[op_ind].counted_pairs.end(); i++) {
				ii++;
			}
		}
		return ii;
	}

	// get a list of all the particle pairs involved in any of the hb order parameters
	vector_of_pairs get_hb_particle_list() {
		vector_of_pairs inds;
		for(int i = 0; i < _hb_parameters_count; i++) {
			for(set_of_pairs::iterator j = _hb_parameters[i].counted_pairs.begin(); j != _hb_parameters[i].counted_pairs.end(); j++) {
				int p_ind = (*j).first;
				int q_ind = (*j).second;
				inds.push_back(std::make_pair(p_ind, q_ind));
			}
		}
		return inds;
	}

	/**
	 * @brief This returns values of ops (indexed by their id_numbers)
	 *
	 * @param param_id
	 * @return values of ops (indexed by ther id_numbers)
	 */
	int get_hb_parameter(int param_id) {
		if(param_id < _hb_parameters_count)
		return _hb_parameters[param_id].get_state();
		else
		return -1;
	}

	double get_hb_cutoff(int param_id) {
		if(param_id < _hb_parameters_count)
		return _hb_parameters[param_id].get_cutoff();
		else
		return HB_CUTOFF;
	}

	double get_distance_parameter(int param_id) {
		if(param_id < _distance_parameters_count) {
			return _distance_parameters[param_id].get_state();
		}
		else {
			return -1;
		}
	}

	/**
	 * @brief For convenience: mapping names to integer ids of the hb order parameters
	 *
	 * @param name
	 * @return
	 */
	int get_hbpar_id_from_name(const char *name);
	const std::string get_name_from_hb_id(int id) {
		return _hb_parameters[id].get_name();
	}

	/**
	 * @brief For convenience: mapping names to integer ids of the distpar order parameters
	 *
	 * @param name
	 * @return
	 */
	int get_distpar_id_from_name(const char *name);
	const std::string get_name_from_distance_id(int id) {
		return _distance_parameters[id].get_name();
	}

	/**
	 * @brief Returns current values of hb_states;
	 * @return
	 *
	 * @warning the returned array can get
	 * @warning overwritten by the subsequent function
	 */
	int *get_hb_states();

	int *get_all_states();

	/// return maximal values that the bond
	/// order parameters can have; warning,
	/// the array can get overwritten by
	/// the previous function
	int *get_max_hb_states();
	int *get_state_sizes();
	double *get_distance_states();

	void print();
	void sprintf_names_and_values(char * str);

	/// to be called in order to evaluate order_parameters:
	/// calculates all minimal distances from list of particles

	void fill_distance_parameters(std::vector<BaseParticle *> &particles, BaseBox * box) {
		for(int i = 0; i < _distance_parameters_count; i++) {
			//_distance_parameters[i].calculate_value(particles, box_side);
			// the next line also updates current_value
			_distance_parameters[i].calculate_state(particles, box);
		}
	}

	/// adds given bonded pair to all values of order parameters
	void add_hb(int a, int b, double energy = -16);
	void remove_hb(int a, int b);

	/// to be called in MC or MD cycles functions
	void reset(); ///sets all hb order params to 0
	void store(); ///stores values of all ops
	void restore(); ///reloads saved values of all ops

	/**
	 * @brief Load ops from file.
	 *
	 * @param _external_filename
	 * @param particles
	 * @param max_N
	 * @return
	 */
	int init_from_file(const char *_external_filename, std::vector<BaseParticle *> &particles, int max_N);

	/**
	 * @brief Sets the log level;
	 * @return
	 */
	void set_log_level (int arg) {_log_level = arg;}

	virtual ~OrderParameters();
};

#endif /* ORDERPARAMETERS_H_ */

