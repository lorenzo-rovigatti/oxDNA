/*
 * HBList.cpp
 *
 *  Created on: Apr 11, 2013
 *      Author: Ben Snodin
 */

#include "HBList.h"

#include <sstream>

template<typename number>
HBList<number>::HBList() {

}

template<typename number>
HBList<number>::~HBList() {

}

template<typename number> void
HBList<number>::init(ConfigInfo<number> &config_info) {
   BaseObservable<number>::init(config_info);
   if (_read_op){
	   _op.init_from_file(_order_parameters_file, this->_config_info.particles, *(this->_config_info.N));
   }
}

template<typename number> void
HBList<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	_read_op = true;
	bool parameters_loaded = false;
	//first try to load parameters from specific op_file key
	if(getInputString(&my_inp, "order_parameters_file", _order_parameters_file, 0) == KEY_NOT_FOUND)
	{
		if(getInputString(&my_inp, "op_file", _order_parameters_file, 0) != KEY_NOT_FOUND)
		{
			parameters_loaded = true;
		}

	}
	else
	{
		parameters_loaded = true;
	}

	//no success, so try to load it from the main simulation input file
	if(parameters_loaded == false)
	{
		if(getInputString(&sim_inp, "order_parameters_file", _order_parameters_file, 0) == KEY_NOT_FOUND)
		{
			if (getInputString(&sim_inp, "op_file", _order_parameters_file, 0) == KEY_NOT_FOUND)
				{
					OX_LOG(Logger::LOG_INFO, "No order parameters file specified in input file; printing indices of any particle pairs that have a hydrogen bond between them");
					_read_op = false;
				}
		}
	}


}

template<typename number>
bool HBList<number>::is_hbond(int p_ind, int q_ind){
	// return true if particle pair p_ind, q_ind have a hydrogen bond between them
	return false;
}

template<typename number>
std::string HBList<number>::get_output_string(llint curr_step) {
	std::stringstream outstr;
	outstr << "# step " << curr_step << "\n";
	// using an order parameters file
	if (_read_op){
		vector_of_pairs inds = _op.get_hb_particle_list();
		for (vector_of_pairs::iterator i = inds.begin(); i != inds.end(); i++) {
			int p_ind = (*i).first;
			int q_ind = (*i).second;
			BaseParticle<number> *p = this->_config_info.particles[p_ind];
			BaseParticle<number> *q = this->_config_info.particles[q_ind];
			number hb_energy = this->_config_info.interaction->pair_interaction_term(DNAInteraction<number>::HYDROGEN_BONDING, p, q);
			if (hb_energy < HB_CUTOFF){
				outstr << p_ind << " " << q_ind << "\n";
			}
		}
	}
	// checking all particle pairs
	else{
		std::vector<ParticlePair<number> > pairs = this->_config_info.lists->get_potential_interactions();
		typename std::vector<ParticlePair<number> >::iterator it;
		for (it = pairs.begin(); it != pairs.end(); it ++ ) {
			BaseParticle<number> * p = (*it).first;
			BaseParticle<number> * q = (*it).second;
			int p_ind = p->get_index();
			int q_ind = q->get_index();
			number hb_energy = this->_config_info.interaction->pair_interaction_term(DNAInteraction<number>::HYDROGEN_BONDING, p, q);
			if (hb_energy < HB_CUTOFF){
				outstr << p_ind << " " << q_ind << "\n";
			}
		}
	}
	return outstr.str();
}

template class HBList<float>;
template class HBList<double>;
