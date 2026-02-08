/*
 * LR_matrix.cpp
 *
 *  Created on: 14 ott 2019
 *      Author: lorenzo
 */

#include "LR_matrix.h"

#include <iostream>

std::ostream &operator<<(std::ostream &os, LR_matrix const m) {
    return os << "[(" << m.v1.x << ", " << m.v1.y << ", " << m.v1.z << "),\n" << 
                    " (" << m.v2.x << ", " << m.v2.y << ", " << m.v2.z << "),\n" <<
                    " (" << m.v3.x << ", " << m.v3.y << ", " << m.v3.z << ")]";
}
