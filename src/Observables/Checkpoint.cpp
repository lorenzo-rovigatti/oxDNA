/*
 * Checkpoint.cpp
 *
 *  Created on: May 16, 2015
 *      Author: flavio
 */

#include <sstream>

#include "Checkpoint.h"

Checkpoint::Checkpoint() {

}

Checkpoint::~Checkpoint() {

}

void Checkpoint::init(ConfigInfo &config_info) {
	_conf.init(config_info);
	this->_config_info = config_info;
}

std::string Checkpoint::get_output_string(llint curr_step) {
	this->_config_info.lists->global_update(true);
	return _conf.get_output_string(curr_step);
}
