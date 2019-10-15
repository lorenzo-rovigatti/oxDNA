/*
 * all_vectors.cpp
 *
 * Created on Apr 26, 2019
 *       Author: Poppleton
 */

#include "AllVectors.h"
#include <sstream>

AllVectors::AllVectors() {

}

AllVectors::~AllVectors() {

}

void AllVectors::init(ConfigInfo &config_info) {
	BaseObservable::init(config_info);
}

std::string AllVectors::get_output_string(llint curr_step) {
	int n = *this->_config_info.N;

	int k = 0;
	std::stringstream outstr;
	for(int i = 0; i < n; i++) {
		for(int j = i + 1; j < n; j++) {
			LR_vector dist;
			LR_vector p1_com, p2_com;
			p1_com = this->_config_info.box->get_abs_pos(this->_config_info.particles[i]);
			p2_com = this->_config_info.box->get_abs_pos(this->_config_info.particles[j]);
			LR_vector vec(0, 0, 0);
			vec = this->_config_info.box->min_image(p1_com, p2_com);
			outstr << vec.x << " " << vec.y << " " << vec.z << "\n";
			k++;
		}
	}
	return outstr.str();
}
