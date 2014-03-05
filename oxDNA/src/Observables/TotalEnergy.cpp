/*
 * TotalEnergy.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#include <sstream>

#include "TotalEnergy.h"

// Initialization of static data members
template<typename number>
llint TotalEnergy<number>::_saved_on = -1;

template<typename number>
number TotalEnergy<number>::_U = (number) 0.;

template<typename number>
number TotalEnergy<number>::_K = (number) 0.;

template<typename number>
TotalEnergy<number>::TotalEnergy() {

}

template<typename number>
TotalEnergy<number>::~TotalEnergy() {

}

template<typename number>
void TotalEnergy<number>::init(ConfigInfo<number> &config_info) {
	_pot_energy.init(config_info);
	_kin_energy.init(config_info);
}

template<typename number>
string TotalEnergy<number>::get_output_string(llint curr_step) {
	number U = get_U(curr_step);
	number K = get_K(curr_step);

	return Utils::sformat("% 10.6lf % 10.6lf % 10.6lf", U, K, U+K);
}

template<typename number>
number TotalEnergy<number>::get_U(llint curr_step) {
	if(curr_step > _saved_on) TotalEnergy::_U = _pot_energy.get_potential_energy();

	return TotalEnergy::_U;
}

template<typename number>
number TotalEnergy<number>::get_K(llint curr_step) {
	if(curr_step > _saved_on) TotalEnergy::_K = _kin_energy.get_kinetic_energy();

	return TotalEnergy::_K;
}

template class TotalEnergy<float>;
template class TotalEnergy<double>;
