/*
 * contact_map.cpp
 *
 * Created on Mar 25, 2019
 *       Author: Poppleton
 */

#include "ContactMap.h"
#include <sstream>

ContactMap::ContactMap() {

}

ContactMap::~ContactMap() {

}

void ContactMap::init(ConfigInfo &config_info) {
	BaseObservable::init(config_info);
}

std::string ContactMap::get_output_string(llint curr_step) {
	int n = *this->_config_info.N;
	int s = ((n * n) - n) / 2;
	number *cmap = new number[s];

	int k = 0;
	for(int i = 0; i < n; i++) {
		for(int j = i + 1; j < n; j++) {
			LR_vector dist;
			LR_vector p1_com, p2_com;
			p1_com = this->_config_info.box->get_abs_pos(this->_config_info.particles[i]);
			p2_com = this->_config_info.box->get_abs_pos(this->_config_info.particles[j]);
			cmap[k] = this->_config_info.box->min_image(p1_com, p2_com).module();
			k++;
		}
	}
	std::stringstream outstr;
	for(int i = 0; i < s; i++) {
		outstr << cmap[i];
		outstr << ' ';
	}

	delete[] cmap;
	return outstr.str();
}
