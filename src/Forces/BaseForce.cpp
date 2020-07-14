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
	_F0 = -1.;
	_rate = -1.;
	_direction = LR_vector(1., 0., 0.);
	_pos0 = LR_vector(0., 0., 0.);
	_site = -1;
	_stiff = 0.;
	_p_ptr = P_VIRTUAL;
}

BaseForce::~BaseForce() {

}
