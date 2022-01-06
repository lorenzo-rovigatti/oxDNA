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
