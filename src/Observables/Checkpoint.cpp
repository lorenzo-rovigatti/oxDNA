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

std::string Checkpoint::get_output_string(llint curr_step) {
	_config_info->lists->global_update(true);
	return _conf.get_output_string(curr_step);
}
