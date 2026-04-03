/*
 * OrderParameters.cpp
 *
 *  Created on: Feb 10, 2012
 *      Author: sulc
 */

#include "OrderParameters.h"

#include <cstdlib>
#include <fstream>

HBParameter::HBParameter() {
	current_value = 0;
	stored_value = 0;
	name = "unnamed_parameter";
	cutoff = HB_CUTOFF;
}

int HBParameter::plus_pair(base_pair& bp,double energy) {
	set_of_pairs::iterator i;
	i = counted_pairs.find(bp);
	if (i != counted_pairs.end()) {
		if(energy < cutoff) // energy is by defualt -16, and if this function was called wihtout speicfying enrgy (default is -16), then it is assumed that the given pair is bonded
		{
			current_value ++;
			return 1;
		}
		else return 0;  //this means that the energy does not qualify for bonding
		
	}
	else return 0;
}

void HBParameter::save_pair_as_parameter(int a, int b) {
	base_pair bp(a, b);
	if (a > b) {
		bp.second = a;
		bp.first = b;
	}
	counted_pairs.insert(bp);
}

int HBParameter::minus_pair(base_pair& bp) {
	set_of_pairs::iterator i;
	i = counted_pairs.find(bp);
	if (i != counted_pairs.end()) {
		current_value --;
		if (current_value < 0) { //this should not happen!
			return -10;
		}
		return -1;
	}
	else return 0;
}

MinDistanceParameter::MinDistanceParameter() {
	state_index = -1;
	stored_state_index = -1;
	n_states = 0;
	current_value = -1.;
	stored_value = -1;
	name = "unnamed_parameter";
}

void MinDistanceParameter::save_pair_as_parameter(int a, int b) {
	base_pair bp(a, b);
	if (a > b) {
		bp.second = a;
		bp.first = b;
	}
	counted_pairs.push_back(bp);
}

int MinDistanceParameter::calculate_state(std::vector<BaseParticle *> &particle_list, BaseBox * box) {
	// if mindistance
	if(_sub_type == 0) {
		LR_vector dist;
		current_value = -1;
		double candidate;
		for(vector_of_pairs::iterator i = counted_pairs.begin(); i != counted_pairs.end(); i++) {
			BaseParticle *p = particle_list[(*i).first];
			BaseParticle *q = particle_list[(*i).second];
			if(p->strand_id == q->strand_id) dist = q->pos - p->pos;
			else dist = box->min_image(p, q);

			if(_use_COM == false) {
				dist += q->int_centers[DNANucleotide::BASE];
				dist -= p->int_centers[DNANucleotide::BASE];
			}
			candidate = dist * dist;

			if(current_value < 0 || candidate < current_value) current_value = candidate;
		}
		current_value = sqrt(current_value);
	}
	// twist
	// Algorithm described in page 26 of Christian Matek's thesis - Statistical mechanics of Nucleic Acids under Mechanical Stress
	else if(_sub_type == 1) {
		double twist = 0;
		for(vector_of_pairs::iterator i = counted_pairs.begin() + 1; i != counted_pairs.end(); i++) {
			// get the position vectors of the current and previous base-pairs
			BaseParticle *first_curr = particle_list[(*i).first];
			BaseParticle *second_curr = particle_list[(*i).second];
			vector_of_pairs::iterator prev = i - 1;
			BaseParticle *first_prev = particle_list[(*prev).first];
			BaseParticle *second_prev = particle_list[(*prev).second];

			// base-pair vectors - vectors that connect the center of mass of the two nucleotides in a base-pair
			LR_vector bp_curr = first_curr->pos - second_curr->pos;
			LR_vector bp_prev = first_prev->pos - second_prev->pos;
			// versor that connects the two base-pair vectors
			LR_vector conn = (bp_curr - bp_prev);
			conn.normalize();
			// project the two base-pair vectors onto the plane normal to the conn vector
			bp_curr -= (conn * bp_curr) * bp_curr;
			bp_prev -= (conn * bp_prev) * bp_prev;
			// finally, compute the angle between the projections and sum it up
			twist += LRACOS(bp_curr * bp_prev);
		}
		//return the twist, which is defined as the number of turns (i.e. the angle in radians divided by 2pi).
		current_value = twist / 2 * M_PI;
	}
	// unknown order parameter.
	else {
		throw oxDNAException("Unknown order parameter subtype %d", _sub_type);
	}

	if(n_states == 0) {
		state_index = -1;
		return state_index;
	}

	if(n_states == 1) {
		state_index = 0;
		return state_index;
	}

	int c = 0;
	while(c < (n_states - 1) && current_value > interfaces[c])
		c++;
	state_index = c;

	return state_index;
}

OrderParameters::OrderParameters() {
	_hb_parameters_count = 0;
	_distance_parameters_count = 0;
	_all_states_count = 0;
	_hb_states = NULL;
	_distance_states = NULL;
	_all_states = NULL;
	_log_level = Logger::LOG_INFO;
}

//adds given bonded pair to all values of order parameters
void OrderParameters::add_hb(int a, int b,double energy) {
	base_pair npair(a, b);
	if (a > b) {
		npair.first = b;
		npair.second = a;
	}
	for (int i = 0; i < _hb_parameters_count; i++) {
		_hb_parameters[i].plus_pair(npair,energy);
	}
}

void OrderParameters::remove_hb(int a, int b) {
	base_pair npair(a, b);
	if (a > b) {
		npair.first = b;
		npair.second = a;
	}
	for (int i = 0; i < _hb_parameters_count; i++) {
		_hb_parameters[i].minus_pair(npair);
	}
}

void OrderParameters::reset(void) {
	for (int i = 0; i < _hb_parameters_count; i++) {
		_hb_parameters[i].reset();
	}
	for (int i = 0; i < _distance_parameters_count; i++) {
		_distance_parameters[i].reset();
	}
}

void OrderParameters::store(void) {
	for (int i = 0; i < _hb_parameters_count; i++) {
		_hb_parameters[i].store_current_value();
	}
	for (int i = 0; i < _distance_parameters_count; i++) {
		_distance_parameters[i].store_current_value();
	}
}

void OrderParameters::restore(void) {
	for (int i = 0; i < _hb_parameters_count; i++) {
		_hb_parameters[i].restore_current_value();
	}
	for (int i = 0; i < _distance_parameters_count; i++) {
		_distance_parameters[i].restore_current_value();
	}
}

int OrderParameters::get_hbpar_id_from_name(const char *name) {
	string stname(name);
	if (_hb_parnames.count(stname) != 0)
		return _hb_parnames[stname];
	else
		return -1;
}

int OrderParameters::get_distpar_id_from_name(const char *name) {
	string stname(name);

	if (_dist_parnames.count(stname) != 0)
		return _dist_parnames[stname];
	else
		return -1;
}

int *OrderParameters::get_state_sizes(void) {
	int i = 0, j = 0;
	
	// the +1 is needed in the following line...
	for (i = 0; i < _hb_parameters_count; i++) _all_states[i] = _hb_parameters[i].get_max_state() + 1;
	
	// but not in the next
	for (j = 0; j < _distance_parameters_count; j++) _all_states[i + j] = _distance_parameters[j].n_states;

	return _all_states;
}


int *OrderParameters::get_max_hb_states(void) {
	for (int i = 0; i < _hb_parameters_count; i++)
		_hb_states[i] = _hb_parameters[i].get_max_state();

	return _hb_states;
}


int *OrderParameters::get_hb_states(void) {
	for (int i = 0; i < _hb_parameters_count; i++)
		_hb_states[i] = _hb_parameters[i].get_state();

	return _hb_states;
}


double * OrderParameters::get_distance_states(void) {
	for (int i = 0; i < _distance_parameters_count; i++)
		_distance_states[i] = _distance_parameters[i].get_state();

	return _distance_states;
}

int *OrderParameters::get_all_states(void) {
	int i = 0, j = 0;
	for (i = 0; i < _hb_parameters_count; i++) _all_states[i] = _hb_parameters[i].get_state();
	for (j = 0; j < _distance_parameters_count; j++) _all_states[i + j] = _distance_parameters[j].get_state_index();
	return _all_states;
}

void OrderParameters::print() {
	int *states;
	states = get_hb_states();
	for(int i = 0; i < _hb_parameters_count; i++) {
		printf("%d ", states[i]);
	}
	double * dists = get_distance_states();
	for(int i = 0; i < _distance_parameters_count; i++)
		printf("%g ", dists[i]);
	printf("\n");
}

void OrderParameters::sprintf_names_and_values(char * str) {
	char tmp[1024];
	*str = '\0';
	int * states = get_hb_states();
	for(int i = 0; i < _hb_parameters_count; i++) {
		sprintf(tmp, "%s: %d; ", get_name_from_hb_id(i).c_str(), states[i]);
		strcat(str, tmp);
	}
	double * dists = get_distance_states();
	for(int i = 0; i < _distance_parameters_count; i++) {
		sprintf(tmp, "%s: %12.9f; ", get_name_from_distance_id(i).c_str(), dists[i]);
		strcat(str, tmp);
	}
	//sprintf("\n");
}

OrderParameters::~OrderParameters() {
	delete[] _hb_states;
	delete[] _distance_states;
	delete[] _all_states;
}

int OrderParameters::init_from_file(const char *_external_filename, std::vector<BaseParticle *> &particles, int max_N) {
	OX_LOG(_log_level, "Parsing order parameter file %s", _external_filename);

	//char line[512], typestr[512];
	int open, justopen, a;
	ifstream external(_external_filename);

	if (!external.good()) throw oxDNAException("Can't read file '%s'", _external_filename);

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
			if (justopen) throw oxDNAException("Syntax error in '%s': nothing between parentheses", _external_filename);
			open--;
		}
		if (open > 1 || open < 0) throw oxDNAException("Syntax error in '%s': parentheses do not match", _external_filename);
		a = external.get();
	}
	external.clear();
	external.seekg(0, ios::beg);

	a = external.get();
	while (external.good()) {
		while (a != '{' && external.good())
		a = external.get();
		if (!external.good()) break;

		a = external.get();
		std::string input_string;
		while (a != '}' && external.good()) {
			input_string += a;
			a = external.get();
		}

		input_file input;
		input.init_from_string(input_string);

		std::string type_str;
		std::string name_str;
		int type;
		int sub_type = 0;
		getInputString(&input, "order_parameter", type_str, 1);
		type = -1;
		if (strcmp(type_str.c_str(), "bond") == 0) type = 0;
		if (strcmp(type_str.c_str(), "mindistance") == 0) {
			type = 1;
			sub_type = 0;
		}
		if (strcmp(type_str.c_str(), "twist") == 0) {
			type = 1;
			sub_type = 1;
		}
		if (type != 0 && type != 1) throw oxDNAException("Parameter type %s not implemented. Aborting", type_str.c_str());
		if (type == 0 || type == 1) {
			getInputString(&input, "name", name_str, 1);
		}
		OX_LOG(_log_level, "Order parameter type %s found, processing", type_str.c_str());

		switch (type) {
			case 0: {
				// HB
				HBParameter newpar;
				_hb_parnames[string(name_str)] = _hb_parameters_count;

				double cutoff;
				if (getInputDouble(&input, "cutoff", &cutoff, 0) == KEY_FOUND) {
					newpar.set_cutoff(cutoff);
					OX_LOG(_log_level, "Using custom cutoff (applies to FFS simulations only) %f ", cutoff);
				}

				vector<string> pair_strings;
				int n_keys = getInputKeys(&input, "pair", &pair_strings, 1);

				for (int l = 0; l < n_keys; l ++) {
					string my_value;
					getInputString (&input, pair_strings[l].c_str(), my_value, 1);

					vector<string> my_values = Utils::split (my_value.c_str(), ',');

					if (my_values.size() != 2) throw oxDNAException("Syntax error in the order parameter file. Found invalid value `%s' for key `%s'. Aborting.", my_value.c_str(), pair_strings[l].c_str());

					// we check if the two strings contain '-'
					bool found1 = my_values[0].find('-') != string::npos;
					bool found2 = my_values[1].find('-') != string::npos;
					if (found1 != found2) throw oxDNAException("Syntax error in the order parameter file. Found invalid value `%s' for key `%s'. Aborting.", my_value.c_str(), pair_strings[l].c_str());

					// handle the case with the list
					if (found1) {
						vector<string> tmps1 = Utils::split(my_values[0].c_str(), '-');
						vector<string> tmps2 = Utils::split(my_values[1].c_str(), '-');

						if (tmps1.size() != 2 || tmps2.size() != 2) throw oxDNAException ("Syntax error in the order parameter file. Found invalid value `%s' for key `%s'. Aborting", my_value.c_str());

						int start1 = atoi(tmps1[0].c_str());
						int end1 = atoi(tmps1[1].c_str());
						int start2 = atoi(tmps2[0].c_str());
						int end2 = atoi(tmps2[1].c_str());
						//printf ("@@ start1 %d, end1 %d, start2 %d, end2 %d\n", start1, end1, start2, end2);

						if (start1 - end1 != start2 - end2 || start1 > end1 || start2 > end2) throw oxDNAException ("Syntax error (d) in the order parameter file. Found invalid list `%s'. Perhaps the indexes are not correct?", my_value.c_str(), pair_strings[l].c_str());

						int cntr = 0;
						for (int m = start1; m <= end1; m ++) {
							int id1 = m;
							int id2 = end2 - cntr;
							if (id1 > max_N || id2 >= max_N) throw oxDNAException("Particle index (%d) out of range while parsing order parameters. Aborting", (id1 > id2) ? id1 : id2);
							if (particles[id1]->btype + particles[id2]->btype != 3) OX_LOG(Logger::LOG_WARNING, "HB pair %d %d not complementary, but still an order parameter", id1, id2);
							newpar.save_pair_as_parameter (id1, id2);
							OX_LOG (_log_level, "--> Adding HB pair (%d, %d) from list `%s' to order parameter `%s'", id1, id2, my_value.c_str(), name_str.c_str());
							cntr ++;
						}
					}
					// handle the normal case, with pairX = <int>, <int>
					else {
						int id1 = atoi (my_values[0].c_str());
						int id2 = atoi (my_values[1].c_str());
						if (id1 > max_N || id2 >= max_N) throw oxDNAException("Particle index (%d) out of range while parsing order parameters. Aborting", (id1 > id2) ? id1 : id2);
						if (particles[id1]->btype + particles[id2]->btype != 3) OX_LOG(Logger::LOG_WARNING, "HB pair %d %d not complementary, but still an order parameter", id1, id2);
						newpar.save_pair_as_parameter (id1, id2);
						OX_LOG (_log_level, "--> Adding HB pair (%d, %d) to order parameter `%s'", id1, id2, name_str.c_str());
					}
				}

				/*
				 int paira, pairb;
				 char strdir[512];
				 char pairname[512];
				 int pairid = 1;
				 sprintf(pairname, "pair%d", pairid);

				 while (getInputString(&input, pairname, strdir, 0) == KEY_FOUND) {
				 //fprintf(stderr, "Loaded %s\n", strdir);
				 tmpi = sscanf(strdir, "%d,%d", &paira, &pairb);

				 if (tmpi != 2)
				 throw oxDNAException("could not parse pairs in HB parameter in order parameters file");
				 if (paira >= max_N || pairb >= max_N)
				 throw oxDNAException("particle index out of range while parsing order parameters");

				 OX_LOG(Logger::LOG_INFO, "--> adding HB pair %d %d ", paira, pairb);

				 newpar.save_pair_as_parameter(paira, pairb);

				 if (particles[paira]->btype + particles[pairb]->btype != 3)
				 OX_LOG(Logger::LOG_WARNING, "HB pair %d %d not complementary, but still an order parameter", paira, pairb);

				 pairid++;
				 sprintf(pairname, "pair%d", pairid);
				 }
				 if (pairid == 1) throw oxDNAException("Error, did not find any parameters to parse. Are pairs correctly numbered?");
				 */

				newpar.set_name (name_str);

				_hb_parameters.push_back(newpar);
				_hb_parameters_count++;

				break;
			}
			case 1: {
				// mindistance
				MinDistanceParameter newpar;
				_dist_parnames[string(name_str)] = _distance_parameters_count;
				newpar._use_COM = false;//defaults to false
				newpar._sub_type = sub_type;

				vector<string> pair_strings;
				int n_keys = getInputKeys(&input, "pair", &pair_strings, 1);

				for (int l = 0; l < n_keys; l ++) {
					string my_value;
					getInputString (&input, pair_strings[l].c_str(), my_value, 1);

					vector<string> my_values = Utils::split (my_value.c_str(), ',');

					if (my_values.size() != 2) throw oxDNAException("Syntax error in the order parameter file. Found invalid value `%s' for key `%s'. Aborting.", my_value.c_str(), pair_strings[l].c_str());

					// we check if the two strings contain '-'
					bool found1 = my_values[0].find('-') != string::npos;
					bool found2 = my_values[1].find('-') != string::npos;
					if (found1 != found2) throw oxDNAException("Syntax error in the order parameter file. Found invalid value `%s' for key `%s'. Aborting.", my_value.c_str(), pair_strings[l].c_str());

					// handle the case with the list
					if (found1) {
						vector<string> tmps1 = Utils::split(my_values[0].c_str(), '-');
						vector<string> tmps2 = Utils::split(my_values[1].c_str(), '-');

						if (tmps1.size() != 2 || tmps2.size() != 2) throw oxDNAException ("Syntax error in the order parameter file. Found invalid value `%s' for key `%s'. Aborting", my_value.c_str());

						int start1 = atoi(tmps1[0].c_str());
						int end1 = atoi(tmps1[1].c_str());
						int start2 = atoi(tmps2[0].c_str());
						int end2 = atoi(tmps2[1].c_str());
						printf ("@@ start1 %d, end1 %d, start2 %d, end2 %d\n", start1, end1, start2, end2);

						if (start1 - end1 != start2 - end2 || start1 > end1 || start2 > end2) throw oxDNAException ("Syntax error (d) in the order parameter file. Found invalid list `%s'. Perhaps the indexes are not correct?", my_value.c_str(), pair_strings[l].c_str());

						int cntr = 0;
						for (int m = start1; m <= end1; m ++) {
							int id1 = m;
							int id2 = end2 - cntr;
							if (id1 > max_N || id2 >= max_N) throw oxDNAException("Particle index (%d) out of range while parsing order parameters. Aborting", (id1 > id2) ? id1 : id2);
							if (particles[id1]->btype + particles[id2]->btype != 3) OX_LOG(Logger::LOG_WARNING, "HB pair %d %d not complementary, but still an order parameter", id1, id2);
							newpar.save_pair_as_parameter (id1, id2);
							OX_LOG (_log_level, "--> Adding HB pair (%d, %d) from list `%s' to order parameter `%s'", id1, id2, my_value.c_str(), name_str.c_str());
							cntr ++;
						}
					}
					// handle the normal case, with pairX = <int>, <int>
					else {
						int id1 = atoi (my_values[0].c_str());
						int id2 = atoi (my_values[1].c_str());
						if (id1 > max_N || id2 >= max_N) throw oxDNAException("Particle index (%d) out of range while parsing order parameters. Aborting", (id1 > id2) ? id1 : id2);
						if (particles[id1]->btype + particles[id2]->btype != 3) OX_LOG(Logger::LOG_WARNING, "HB pair %d %d not complementary, but still an order parameter", id1, id2);
						newpar.save_pair_as_parameter (id1, id2);
						OX_LOG (_log_level, "--> Adding HB pair (%d, %d) to order parameter `%s'", id1, id2, name_str.c_str());
					}
				}

				/*
				 // parsing of the interfaces...
				 newpar.n_states = 0;
				 if (getInputString (&input, "interfaces", strdir, 0) == KEY_FOUND) {
				 // we need to parse the interfaces; we expect at least one float number.
				 char * aux = strtok (strdir, ",");
				 int i = 0;
				 while (aux != NULL) {
				 double tmpf;
				 sscanf(aux, "%lf", &tmpf);
				 newpar.interfaces.push_back(tmpf);
				 OX_LOG(_log_level, " ---> found interface %i; %lf", i, newpar.interfaces[i]);
				 i ++;
				 aux = strtok (NULL, ",");
				 }
				 newpar.n_states = i + 1; // the number of states is the number of interfaces + 1;
				 }
				 */
				string strdir;
				if (getInputString (&input, "interfaces", strdir, 0) == KEY_FOUND) {
					std::vector<string> interfaces = Utils::split (strdir.c_str(), ',');
					for (unsigned int i = 0; i < interfaces.size(); i ++) {
						newpar.interfaces.push_back (atof(interfaces[i].c_str()));
						OX_LOG(_log_level, " ---> found interface %i; %lf", i, newpar.interfaces[i]);
					}
					newpar.n_states = interfaces.size() + 1;
				}

				//optional use of COM-COM distance instead of base-base distance
				/*
				 char tmpstr2[256];
				 if (getInputString(&input, "use_COM", tmpstr2, 0) == KEY_FOUND) {
				 unsigned int COM=0;
				 sscanf(tmpstr2, "%u", &COM);
				 if(COM==0) { newpar._use_COM = false; }
				 else if(COM==1) {
				 newpar._use_COM = true;
				 OX_LOG (_log_level, "using COM-COM distances instead of base-base");
				 }
				 else { OX_LOG (_log_level, "unsupported value for use_COM"); }
				 }*/

				// optional use of COM-COM distance instead of base-base distance
				if (getInputBool(&input, "use_COM", &newpar._use_COM, 0) == KEY_FOUND) {
					if (newpar._use_COM) OX_LOG (_log_level, "using COM-COM distances instead of base-base");
				}

				newpar.set_name (name_str);
				_distance_parameters.push_back(newpar);
				_distance_parameters_count ++;

				break;

			}

			default:
			OX_LOG(Logger::LOG_WARNING, "Probably should't reach this point. Hoping for the best");
			break;
		}
	}

	if (_hb_parameters_count > 0) _hb_states = new int[_hb_parameters_count];
	if (_distance_parameters_count > 0) _distance_states = new double[_distance_parameters_count];
	_all_states_count = _hb_parameters_count + _distance_parameters_count;
	if (_all_states_count > 0) _all_states = new int[_all_states_count];

	OX_LOG(_log_level, " File %s parsed; found %d hb_dim, %d dist_dim", _external_filename, _hb_parameters_count, _distance_parameters_count);

	return 0;
}

