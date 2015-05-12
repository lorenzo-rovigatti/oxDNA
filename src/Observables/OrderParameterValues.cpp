/*
 * OrderParameterValues.cpp
 *
 *  Created on: Mar 2, 2013
 *      Author: petr
 */

#include "OrderParameterValues.h"

#include <sstream>

template<typename number>
OrderParameterValues<number>::OrderParameterValues() {

}

template<typename number>
OrderParameterValues<number>::~OrderParameterValues() {

}

template<typename number> void
OrderParameterValues<number>::init(ConfigInfo<number> &config_info) {
   BaseObservable<number>::init(config_info);

   ifstream tmpf(_order_parameters_file);
   if(!tmpf.good ()) throw oxDNAException ("(OrderParameterValues.cpp) Can't read file '%s'", _order_parameters_file);

   _op.init_from_file(_order_parameters_file, this->_config_info.particles, *(this->_config_info.N));

}

template<typename number> void
OrderParameterValues<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	bool parameters_loaded = false;
	//first try load parameters from specific op_file key
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

	//no success, so try it load from the main simulation input file
	if(parameters_loaded == false)
	{
		if(getInputString(&sim_inp, "order_parameters_file", _order_parameters_file, 0) == KEY_NOT_FOUND)
		{
			getInputString(&sim_inp, "op_file", _order_parameters_file, 1);
		}
	}


}

template<typename number>
std::string OrderParameterValues<number>::get_output_string(llint curr_step) {

	_op.reset();
	_op.fill_distance_parameters(this->_config_info.particles,*this->_config_info.box_side);
	
	std::vector<ParticlePair<number> > neighbour_pairs = this->_config_info.lists->get_potential_interactions();
	
	for (int i = 0; i < (int)neighbour_pairs.size(); i++) {
		BaseParticle<number> * p = neighbour_pairs[i].first;
		BaseParticle<number> * q = neighbour_pairs[i].second;

		number hb_energy = this->_config_info.interaction->pair_interaction_term(DNAInteraction<number>::HYDROGEN_BONDING, p, q);
		if(hb_energy < HB_CUTOFF) _op.add_hb(p->index,q->index);
	}

	std::stringstream outstr;
	for (int i = 0; i < _op.get_hb_parameters_count(); i++)
	{
		outstr << _op.get_hb_parameter(i) << " " ;
	}
	for (int i = 0; i < _op.get_distance_parameters_count(); i++)
	{
		outstr << _op.get_distance_parameter(i) << " " ;
	}
	//outstr << '\n';

	return outstr.str();
}

template class OrderParameterValues<float>;
template class OrderParameterValues<double>;
