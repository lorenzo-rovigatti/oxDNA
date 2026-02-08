/*
 * LR_vector.cpp
 *
 *  Created on: 14 ott 2019
 *      Author: lorenzo
 */

#include "LR_vector.h"

#include <iostream>

LR_vector operator*(const number S, const LR_vector &v) {
	return v * S;
}

std::ostream &operator<<(std::ostream &stream, const LR_vector &vector) {
	stream << "(" << vector.x << ", " << vector.y << ", " << vector.z << ")";
	return stream;
}