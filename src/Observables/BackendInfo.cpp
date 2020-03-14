/*
 * BackendInfo.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#include "BackendInfo.h"

BackendInfo::BackendInfo() {

}

BackendInfo::~BackendInfo() {

}

std::string BackendInfo::get_output_string(llint curr_step) {
	return *(_config_info->backend_info);
}
