/*
 * Density.cpp
 *
 *  Created on: Oct 30, 2013
 *      Author: Flavio
 */

#include "Density.h"
#include "../Utilities/Utils.h"

template<typename number>
Density<number>::Density() {

}

template<typename number>
Density<number>::~Density() {

}

template<typename number>
std::string Density<number>::get_output_string(llint curr_step) {
	int N = *this->_config_info.N;
	number V = pow (*this->_config_info.box_side, 3);
	return Utils::sformat ("%# 10.8g", N / V);
}

template<typename number>
void Density<number>::get_settings (input_file &my_inp, input_file &sim_inp) {

}

template class Density<float>;
template class Density<double>;
