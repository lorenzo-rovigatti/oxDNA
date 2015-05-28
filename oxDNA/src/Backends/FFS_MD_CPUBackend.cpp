/*
 * MD_CPUBackend.cpp
 *
 *  Created on: 03/set/2010
 *      Author: sulc
 */

//#include "SimManager.h"
#include "FFS_MD_CPUBackend.h"
#include <sstream>
#include <cctype>

#include "../Interactions/DNAInteraction.h"
#include "../Interactions/DNA2Interaction.h"
#include "../Interactions/RNAInteraction2.h"
#include "../Managers/SimManager.h"

using namespace std;

bool parsed_expression::parse_expression(const char *expression, OrderParameters *op) {

	string s(expression);
	//s.erase(remove_if(s.begin(), s.end(),' '), s.end()); //remove spaces
	s.erase(remove_if(s.begin(), s.end(), static_cast<int(*)(int)>( isspace )), s.end());

	//cout << "Processing " << s << endl;
	int condition_type = -1;
	int pos = 0;
	unsigned int i;
	for (i = 0; i < s.length() - 1; i++) {
		if (s[i] == '<') {
			if (s[i + 1] == '=') {
				condition_type = LEQ;
				pos = i + 2;
				break;
			} else {
				condition_type = LT;
				pos = i + 1;
				break;
			}
		} else if (s[i] == '>') {
			if (s[i + 1] == '=') {
				condition_type = GEQ;
				pos = i + 2;
				break;
			} else {
				condition_type = GT;
				pos = i + 1;
				break;
			}

		} else if (s[i] == '=' && s[i + 1] == '=') {
			condition_type = EQ;
			pos = i + 2;
			break;
		}

	}
	if (condition_type == -1) {
		return false;
	}
	compare_type = condition_type;

	//cout << "Processed as " << condition_type << " and left with  "
	//		<< s.substr(0, i) << endl;

	string parname = s.substr(0, i);
	parameter_index = op->get_hbpar_id_from_name(parname.c_str());
	expression_type = 0;
	if (parameter_index == -1) { //it is a distance parameter
		parameter_index = op->get_distpar_id_from_name(parname.c_str());
		expression_type = 1;
	}
	if (parameter_index == -1) {
		//cout << "did not find parameter " << parname.c_str() << endl;
		return false;
	}

	stringstream str;
	str << s.substr(pos);

	if (expression_type == 0)
		str >> value_to_compare_with_hb;
	else if (expression_type == 1)
		str >> value_to_compare_with_dist;

	return true;

}

bool parsed_expression::eval_expression(OrderParameters *op) {

	if (expression_type == 0) {

		int current_value = op->get_hb_parameter(parameter_index);
		//cout << "Evaluating HB with " << current_value << " to be " << expression_type << " with " << value_to_compare_with_hb << endl;
		switch (compare_type) {
		case GEQ:
			return current_value >= value_to_compare_with_hb ? true : false;
			break;

		case LEQ:
			return current_value <= value_to_compare_with_hb ? true : false;
			break;

		case LT:
			return current_value < value_to_compare_with_hb ? true : false;
			break;
		case GT:
			return current_value > value_to_compare_with_hb ? true : false;
			break;
		case EQ:
			return current_value == value_to_compare_with_hb ? true : false;
			break;
		default:
			break;
		}
	} else if (expression_type == 1) {
		double current_value = op->get_distance_parameter(parameter_index);
		switch (compare_type) {
		case GEQ:
			return current_value >= value_to_compare_with_dist ? true : false;
			break;

		case LEQ:
			return current_value <= value_to_compare_with_dist ? true : false;
			break;

		case LT:
			return current_value < value_to_compare_with_dist ? true : false;
			break;
		case GT:
			return current_value > value_to_compare_with_dist ? true : false;
			break;
		case EQ:
			return current_value == value_to_compare_with_dist ? true : false;
			break;
		default:
			break;
		}
	} else {
		return false;
	}

	return false;

}

bool parsed_condition::parse_condition(const char *expression, OrderParameters *op)
{
	//input_file *obs_input = Utils::get_input_file_from_string(obs_string);
	string str(expression);

	//remove { and }
	str.erase(std::remove(str.begin(), str.end(), '{'), str.end());
	str.erase(std::remove(str.begin(), str.end(), '}'), str.end());

	//remove spaces
	str.erase(remove_if(str.begin(), str.end(), static_cast<int(*)(int)>( isblank )), str.end());



	std::stringstream input(str.c_str());
	std::string str_expression;

	while(std::getline(input,str_expression,'\n'))
	{
	    parsed_expression pe;

	    if(str_expression.length() > 1)
	    {
	      //printf("Parsing %s\n",str_expression.c_str());
	      if (! pe.parse_expression(str_expression.c_str(),  op ))
	      {
	    	 return false;
	      }
	      else
	      {
	     	all_expressions.push_back(pe);
	      }
	    }

	}


	return true;
}

bool parsed_condition::eval_condition(OrderParameters *op)
{
	for (vector<parsed_expression>::iterator i = this->all_expressions.begin();
					i != this->all_expressions.end(); i++) {
				if (!(*i).eval_expression(op) )
					return false;
	}

	return true;

}

template<typename number>
FFS_MD_CPUBackend<number>::FFS_MD_CPUBackend() : MD_CPUBackend<number>() {
	this->_is_CUDA_sim = false;
}

template<typename number>
FFS_MD_CPUBackend<number>::~FFS_MD_CPUBackend() {

}

template<typename number>
void FFS_MD_CPUBackend<number>::get_settings(input_file &inp) {
	MD_CPUBackend<number>::get_settings(inp);
	getInputString(&inp, "order_parameters_file", _order_parameters_file, 1);
	getInputString(&inp, "ffs_file", _ffs_file, 1);
}

template<typename number>
number FFS_MD_CPUBackend<number>::pair_interaction_nonbonded_DNA_with_op(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	if(r->norm() >= this->_sqr_rcut) return (number) 0.f;

	number energy = this->_interaction->pair_interaction_term(DNAInteraction<number>::HYDROGEN_BONDING,p, q, r, update_forces);
	if (energy <= MAX_BOND_CUTOFF) {
		this->_op.add_hb(q->index, p->index,energy);
	}
	energy += this->_interaction->pair_interaction_term(DNAInteraction<number>::NONBONDED_EXCLUDED_VOLUME,p, q, r, update_forces);
	energy += this->_interaction->pair_interaction_term(DNAInteraction<number>::CROSS_STACKING,p, q, r, update_forces);

	// all interactions except DNA2Interaction use the DNAInteraction coaxial stacking
        if(dynamic_cast<DNA2Interaction<number> *>(this->_interaction) == NULL) energy += this->_interaction->pair_interaction_term(DNAInteraction<number>::COAXIAL_STACKING, p, q, r, update_forces);
	
        if(dynamic_cast<DNA2Interaction<number> *>(this->_interaction) != NULL) {
		energy += this->_interaction->pair_interaction_term(DNA2Interaction<number>::COAXIAL_STACKING, p, q, r, update_forces);
		energy += this->_interaction->pair_interaction_term(DNA2Interaction<number>::DEBYE_HUCKEL, p, q, r, update_forces);
	}
	else if(dynamic_cast<RNA2Interaction<number> *>(this->_interaction) != NULL) {
	        energy += this->_interaction->pair_interaction_term(RNA2Interaction<number>::DEBYE_HUCKEL, p, q, r, update_forces);
	}

	return energy;
}

template<typename number>
void FFS_MD_CPUBackend<number>::_ffs_compute_forces(void) {
	this->_U = this->_U_hydr = (number) 0;
	for (int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = this->_particles[i];
		this->_U += this->_interaction->pair_interaction_bonded(p, P_VIRTUAL, NULL, true);

		std::vector<BaseParticle<number> *> neighs = this->_lists->get_neigh_list(p);
		for(unsigned int n = 0; n < neighs.size(); n++) {
			BaseParticle<number> *q = neighs[n];
			// if implemented implicitely in the lists
			//if (p->index < q->index) this->_U += this->pair_interaction_nonbonded_DNA_with_op(p, q, NULL, true);
			this->_U += this->pair_interaction_nonbonded_DNA_with_op(p, q, NULL, true);
		}
	}
}

template<typename number>
void FFS_MD_CPUBackend<number>::sim_step(llint curr_step) {
	this->_mytimer->resume();

	this->_op.reset();

	this->_timer_first_step->resume();
	this->_first_step(curr_step);
	this->_timer_first_step->pause();

	this->_timer_lists->resume();
	if(!this->_lists->is_updated()) {
		this->_lists->global_update();
		this->_N_updates++;
	}
	this->_timer_lists->pause();

	this->_timer_forces->resume();
	this->_ffs_compute_forces();
	this->_second_step();
	this->_timer_forces->pause();

	this->_timer_thermostat->resume();
	this->_thermostat->apply (this->_particles, curr_step);
	this->_timer_thermostat->pause();

	_op.fill_distance_parameters<number>(this->_particles, this->_box_side);

	//cout << "I just stepped and bond parameter is " << _op.get_hb_parameter(0) << " and distance is " << _op.get_distance_parameter(0) << endl;
	if (this->check_stop_conditions()) {
		SimManager::stop = true;
		OX_LOG(Logger::LOG_INFO, "Reached stop conditions, stopping in step %lld", curr_step);
		char tmp[1024];
		_op.sprintf_names_and_values(tmp);
		OX_LOG(Logger::LOG_INFO, "FFS final values: %s", tmp);
	}

	this->_mytimer->pause();
}


template<typename number>
void FFS_MD_CPUBackend<number>::init() {
	MD_CPUBackend<number>::init();

	this->_sqr_rcut = this->_interaction->get_rcut() * this->_interaction->get_rcut();

	_op.init_from_file(_order_parameters_file.c_str(), this->_particles, this->_N);
	init_ffs_from_file(_ffs_file.c_str());
	OX_LOG(Logger::LOG_INFO, "Setting initial value for the order parameter...");
	_op.fill_distance_parameters<number>(this->_particles, this->_box_side);

}

/*
 File format:
 {
 action = stop_and
 condition1 = params_a > 1
 condition2 = params_b >= 4
 condition3 = params_c < 4
 }
 */
template<typename number>
void FFS_MD_CPUBackend<number>::init_ffs_from_file(const char *fname) {

	FILE * fin = 0;
	fin = fopen(fname, "r");
	if (fin == NULL)
		throw oxDNAException("Cannot open %s", fname);

	input_file input;

	loadInput(&input, fin);

	char pairname[512];
	string strexpr;
	int pairid = 1;
	sprintf(pairname, "condition%d", pairid);
	while (getInputString(&input, pairname, strexpr, 0) == KEY_FOUND) {

		parsed_condition newcondition;
		//cout << "Parsing " << strexpr << endl;
		if (!newcondition.parse_condition(strexpr.c_str(), &_op))
			throw oxDNAException("Failed to parse condition %s, please check the file format and parameter names ", pairname);

		this->_conditions.push_back(newcondition);
		pairid++;
		sprintf(pairname, "condition%d", pairid);
	}

	fclose(fin);
}

template<typename number>
char * FFS_MD_CPUBackend<number>::get_op_state_str(void) {
	int * state = _op.get_hb_states();
	char * aux;
	aux = (char *) _state_str;
	for (int i = 0; i < _op.get_hb_parameters_count(); i++) {
		sprintf(aux, "%2d ", state[i]);
		aux = (char *) _state_str + strlen(_state_str);
	}

	double * dstate = _op.get_distance_states();
	for (int i = 0; i < _op.get_distance_parameters_count(); i++) {
		sprintf(aux, "%5.2f ", dstate[i]);
		aux = (char *) _state_str + strlen(_state_str);
	}
	return _state_str;
}

template<typename number>
void FFS_MD_CPUBackend<number>::print_observables(llint curr_step) {
	this->_backend_info = get_op_state_str();
	MDBackend<number>::print_observables(curr_step);
}

template<typename number>
bool FFS_MD_CPUBackend<number>::check_stop_conditions(void) {
 //check if any of the stop conditions is true
		for (vector<parsed_condition>::iterator i = this->_conditions.begin();
				i != this->_conditions.end(); i++) {
			if ((*i).eval_condition(&(this->_op)))
			{
				OX_LOG(Logger::LOG_INFO,"Reached condition%d",1+std::distance(_conditions.begin(), i));
				return true;
			}
		}
		return false;
}




template class FFS_MD_CPUBackend<float> ;
template class FFS_MD_CPUBackend<double> ;
