/*
 * LR_vector.cpp
 *
 *  Created on: 14 ott 2019
 *      Author: lorenzo
 */

#include "LR_vector.h"

LR_vector operator*(const number S, const LR_vector &v) {
	return v * S;
}
