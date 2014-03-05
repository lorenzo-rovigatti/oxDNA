/*
 * BackendInfo.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#include "BackendInfo.h"

template<typename number>
BackendInfo<number>::BackendInfo() {

}

template<typename number>
BackendInfo<number>::~BackendInfo() {

}

template<typename number>
std::string BackendInfo<number>::get_output_string(llint curr_step) {
	return *(this->_config_info.backend_info);
}

template class BackendInfo<float>;
template class BackendInfo<double>;

