/*
 * OrderParameters.h
 *
 *  Created on: Feb 10, 2012
 *      Author: sulc
 */

#ifndef ORDERPARAMETERS_H_
#define ORDERPARAMETERS_H_

#ifndef HB_CUTOFF
#define HB_CUTOFF ((number)(-0.1f))
#endif

#include <set>
#include <utility>
#include <vector>
#include <string>
#include <map>
#include "IOManager.h"
#include "Particle.h"

//using namespace std;

typedef std::pair<int, int> base_pair;
typedef std::vector<base_pair> vector_of_pairs;

struct classcomp {
	bool operator()(const base_pair& a, const base_pair& b) const {
		if (a.first < b.first)
			return true;
		else if (a.first == b.first && a.second < b.second)
			return true;
		else
			return false;
	}
};

typedef std::set<base_pair, classcomp> set_of_pairs;

struct HBParameter {
	set_of_pairs counted_pairs; //this is what defines order parameter: it is loaded from external init file
	int current_value;
	int stored_value;
	std::string name;

	HBParameter(void);

	void save_pair_as_parameter(int a, int b);

	void set_current_value(int new_value) {
		current_value = new_value;
	}
	int store_current_value(void) {
		stored_value = current_value;
		return current_value;
	}
	int restore_current_value(void) {
		current_value = stored_value;
		return current_value;
	}

	int get_state(void) {
		return current_value;
	}

	int get_max_state(void) {
		return counted_pairs.size();
	}

	std::string get_name(void) {
		return name;
	}
	void set_name(std::string newname) {
		name = newname;
	}

	int plus_pair(base_pair& bp);
	int minus_pair(base_pair& bp);
	void reset(void) {
		current_value = 0;
	}

};

//op that saves minimum distance of base distances between pairs
struct MinDistanceParameter {
	vector_of_pairs counted_pairs; //this is what defines order parameter: it is loaded from external init file
	double current_value;
	double stored_value;
	std::string name;

	MinDistanceParameter(void);
	void save_pair_as_parameter(int a, int b);

	void set_current_value(int new_value) {
		current_value = new_value;
	}
	double store_current_value(void) {
		stored_value = current_value;
		return current_value;
	}
	double restore_current_value(void) {
		current_value = stored_value;
		return current_value;
	}

	double get_state(void) {
		return current_value;
	}

	std::string get_name(void) {
		return name;
	}
	void set_name(std::string newname) {
		name = newname;
	}

	template<typename number> double calculate_value(
			Particle<number> *particle_list, double box_side) {
		LR_vector<number> dist;
		current_value = -1;
		double candidate;
		for (vector_of_pairs::iterator i = counted_pairs.begin();
				i != counted_pairs.end(); i++) {
			dist = (particle_list[(*i).first].pos
					+ particle_list[(*i).first].pos_base);
			candidate = (double) dist.sqr_min_image_distance(
					(particle_list[(*i).second].pos
							+ particle_list[(*i).second].pos_base), box_side);
			//candidate = static_cast<double>(dist.norm());
			if (current_value < 0 || candidate < current_value) {
				current_value = candidate;
			}

		}

		current_value = sqrt(current_value);
		return current_value;
	}

};

/*
 Example of order_parameter file:
 {
 order_parameter = bond
 name = all_bonds_01
 pair1 = 14, 20
 pair2 = 14, 16
 pair3 = 17, 41
 }
 {
 order_parameter = mindistance
 name = dist_strand_12
 pair1 = 14,58
 pair2 = 12,45
 }

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

public:
	OrderParameters();

	//access functions
	int get_hb_parameters_count() {
		return _hb_parameters_count;
	}
	int get_distance_parameters_count() {
		return _distance_parameters_count;
	}

	//this returns values of ops (indexed by theri id_numbers)
	int get_hb_parameter(int param_id) {
		if (param_id < _hb_parameters_count)
			return _hb_parameters[param_id].get_state();
		else
			return -1;
	}
	double get_distance_parameter(int param_id) {
		if (param_id < _distance_parameters_count)
			return _distance_parameters[param_id].get_state();
		else
			return -1;
	}

	//for convenience: mapping names to integer ids of the order parameters
	int get_hbpar_id_from_name(const char *name);
	const char *get_name_from_hb_id(int id) {
		return _hb_parameters[id].get_name().c_str();
	}
	int get_distpar_id_from_name(const char *name);
	const char *get_name_from_distance_id(int id) {
		return _distance_parameters[id].get_name().c_str();
	}

	int * get_hb_states(void); //returns current values of hb_states;
							   //warning, the returnd array can get
							   //overwritten by the subsequent function
	int * get_max_hb_states(void); // return maximal values that the bond
								   // order parameters can have; warning,
								   // the arraty can get overwritten by 
								   // the previous function
	double *get_distance_states(void);

	void print(void) {
		int * states;
		states = get_hb_states();
		for (int i = 0; i < _hb_parameters_count; i++) {
			printf("%d ", states[i]);
		}
		printf("\n");
	}

	//to be called in order to evaluate order_parameters:
	//calculates all minimal distances from list of particles
	template<typename number> void fill_distance_parameters(
			Particle<number> *particles, double box_side) {

		for (int i = 0; i < _distance_parameters_count; i++) {
			_distance_parameters[i].calculate_value(particles, box_side);
		}
	}

	//adds given bonded pair to all values of order parameters
	void add_hb(int a, int b);
	void remove_hb(int a, int b);

	//to be called in MC or MD cycles functions
	void reset(void); //sets all hb order params to 0
	void store(void); //stores values of all ops
	void restore(void); // reloads saved values of all ops

	// load ops from file; This unfortunately needs to 
	// be defined in the header file...
	template<typename number>
	int init_from_file(const char *_external_filename,
			Particle<number> * particles, int max_N, IOManager *_IO) {

		_IO->log(_IO->LOG_INFO, "Parsing order parameter file %s",
				_external_filename);

		//char line[512], typestr[512];
		int open, justopen, a;
		ifstream external(_external_filename);

		if (!external.good())
			_IO->die("Can't read file '%s'", _external_filename);

		justopen = open = 0;
		a = external.get();
		int total_parameters = 0;
		while (external.good()) {
			justopen = 0;
			if (a == '{') {
				open++;
				justopen = 1;
				total_parameters++;
			}
			if (a == '}') {
				if (justopen)
					_IO->die(
							"Syntax error in '%s': nothing between parentheses",
							_external_filename);
				open--;
			}
			if (open > 1 || open < 0)
				_IO->die("Syntax error in '%s': parentheses do not match",
						_external_filename);
			a = external.get();
		}
		external.clear();
		external.seekg(0, ios::beg);

		a = external.get();
		while (external.good()) {
			while (a != '{' && external.good())
				a = external.get();
			if (!external.good())
				break;
			// this function create a temporary file which is destroyed
			// upon calling fclose
			// the temporary file is opened with "wb+" flags
			FILE *temp = tmpfile();
			// FILE *temp = fopen(".merda", "w");
			_IO->log(_IO->LOG_INFO, "   Using temporary file");
			a = external.get();
			while (a != '}' && external.good()) {
				fprintf(temp, "%c", a);
				a = external.get();
			}

			rewind(temp);
			input_file input;
			char type_str[512];
			char name_str[512];
			int type;
			loadInput(&input, temp);
			getInputString(&input, "order_parameter", type_str, 1);
			type = -1;
			if (strcmp(type_str, "bond") == 0)
				type = 0;
			if (strcmp(type_str, "mindistance") == 0)
				type = 1;

			if (type != 0 && type != 1)
				_IO->die("Parameter type %s not implemented. Aborting",
						type_str);
			int tmpi;
			if (type == 0 || type == 1) {
				getInputString(&input, "name", name_str, 1);

			}
			_IO->log(_IO->LOG_INFO, "Order parameter type %s found, processing",
					type_str);
			switch (type) {
			case 0: {
				// HB
				HBParameter newpar;
				_hb_parnames[string(name_str)] = _hb_parameters_count;

				int paira, pairb;
				char strdir[512];
				char pairname[512];
				int pairid = 1;
				sprintf(pairname, "pair%d", pairid);
				while (getInputString(&input, pairname, strdir, 0) == KEY_FOUND) {
					//fprintf(stderr, "Loaded %s\n", strdir);
					tmpi = sscanf(strdir, "%d,%d", &paira, &pairb);

					if (tmpi != 2)
						_IO->die(
								"could not parse pairs in HB parameter in order parameters file");
					if (paira >= max_N || pairb >= max_N)
						_IO->die(
								"particle index out of range while parsing order parameters");
					_IO->log(_IO->LOG_INFO, "--> adding HB pair %d %d ", paira,
							pairb);

					newpar.save_pair_as_parameter(paira, pairb);

					if (particles[paira].btype + particles[pairb].btype != 3)
						_IO->log(
								_IO->LOG_WARNING,
								"HB pair %d %d not complementary, but still an order parameter",
								paira, pairb);

					pairid++;
					sprintf(pairname, "pair%d", pairid);
				}
				if (pairid == 1)
					_IO->die(
							"Error, did not find any parameters to parse. Are pairs correctly numbered?");

				_hb_parameters.push_back(newpar);
				_hb_parameters_count++;

				break;
			}
			case 1: {
				// mindistance
				MinDistanceParameter newpar;
				_dist_parnames[string(name_str)] = _distance_parameters_count;

				int paira, pairb;
				char strdir[512];
				char pairname[512];
				int pairid = 1;
				sprintf(pairname, "pair%d", pairid);
				while (getInputString(&input, pairname, strdir, 0) == KEY_FOUND)

				{
					tmpi = sscanf(strdir, "%d,%d", &paira, &pairb);

					if (tmpi != 2)
						_IO->die("could not parse pairs order parameters file");
					if (paira >= max_N || pairb >= max_N)
						_IO->die(
								"particle index out of range while parsing order parameters");
					_IO->log(_IO->LOG_INFO,
							"--> adding mindistance pair %d %d ", paira, pairb);

					newpar.save_pair_as_parameter(paira, pairb);

					if (particles[paira].btype + particles[pairb].btype != 3)
						_IO->log(
								_IO->LOG_WARNING,
								" pair %d %d not complementary, but still an order parameter for minimal distance",
								paira, pairb);
					pairid++;
					sprintf(pairname, "pair%d", pairid);
				}
				if (pairid == 1)
					_IO->die(
							"Error, did not find any parameters to parse. Are pairs correctly numbered?");

				_distance_parameters.push_back(newpar);
				_distance_parameters_count++;
				break;

			}

			default:
				_IO->log(
						_IO->LOG_INFO,
						"Probably should't reach this point. Hoping for the best");
				break;
			}
			cleanInputFile(&input);
			fclose(temp);
		}
		_IO->log(_IO->LOG_INFO, " File %s parsed", _external_filename);

		if (_hb_parameters_count > 0)
			_hb_states = new int[_hb_parameters_count];
		if (_distance_parameters_count > 0)
			_distance_states = new double[_distance_parameters_count];

		return 0;
	}
	virtual ~OrderParameters();
};

#endif /* ORDERPARAMETERS_H_ */

