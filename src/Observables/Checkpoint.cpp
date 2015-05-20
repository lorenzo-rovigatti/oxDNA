/*
 * Checkpoint.cpp
 *
 *  Created on: May 16, 2015
 *      Author: flavio
 */

#include <sstream>

#include "Checkpoint.h"

template<typename number>
Checkpoint<number>::Checkpoint() {

}

template<typename number>
Checkpoint<number>::~Checkpoint() {

}

template<typename number>
void Checkpoint<number>::init(ConfigInfo<number> &config_info) {
	_conf.init(config_info);
	this->_config_info = config_info; 
}

template<typename number>
string Checkpoint<number>::get_output_string(llint curr_step) {
	this->_config_info.lists->global_update(true);
	return _conf.get_output_string(curr_step);
}

template class Checkpoint<float>;
template class Checkpoint<double>;
