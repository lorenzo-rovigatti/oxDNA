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

void ContactMap::init() {
	BaseObservable::init();
}

std::string ContactMap::get_output_string(llint curr_step) {
	int N = _config_info->N();
	int s = ((N * N) - N) / 2;
	static std::vector<number> cmap(s);

	int k = 0;
	for(int i = 0; i < N; i++) {
		for(int j = i + 1; j < N; j++) {
			LR_vector dist;
			LR_vector p1_com, p2_com;
			p1_com = _config_info->box->get_abs_pos(_config_info->particles()[i]);
			p2_com = _config_info->box->get_abs_pos(_config_info->particles()[j]);
			cmap[k] = _config_info->box->min_image(p1_com, p2_com).module();
			k++;
		}
	}
	std::stringstream outstr;
	for(int i = 0; i < s; i++) {
		outstr << cmap[i];
		outstr << ' ';
	}

	return outstr.str();
}
