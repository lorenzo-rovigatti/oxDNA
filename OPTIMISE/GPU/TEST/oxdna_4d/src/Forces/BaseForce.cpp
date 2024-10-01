/*
 * BaseForce.cpp
 *
 *  Created on: 22 giu 2017
 *      Author: lorenzo
 */

#include "BaseForce.h"

#include "../Boxes/BaseBox.h"

#include <vector>
#include <string>
#include <tuple>

BaseForce::BaseForce() {
	_stiff = 0.;
	_p_ptr = P_VIRTUAL;
}

BaseForce::~BaseForce() {

}

std::tuple<std::vector<int>, std::string> BaseForce::init(input_file &inp) {
	getInputString(&inp, "group_name", _group_name, 0);
	getInputString(&inp, "id", _id, 0);
	getInputString(&inp, "type", _type, 1);

	return std::make_tuple(std::vector<int>(), "BaseForce");
}
