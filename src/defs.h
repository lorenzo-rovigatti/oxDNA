/**
 * @file    defs.h
 * @date    13/ott/2009
 * @author  lorenzo
 *
 *
 */

#ifndef DEFS_H_
#define DEFS_H_

#define MAX_EXT_FORCES 15

#define PI 3.141592653589793238462643f
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define LRACOS(x) (((x) > 1) ? (number) 0 : ((x) < -1) ? (number) PI : acos(x))

#define CHECK_BOX(my_class, inp) 	std::string box_type ("");\
	if (getInputString(&inp, "box_type", box_type, 0) == KEY_FOUND) {\
		if (box_type.compare("cubic") != 0) \
			throw oxDNAException ("%s only works with cubic box! Aborting", my_class);\
	}\

#define P_A 0
#define P_B 1
#define P_VIRTUAL (NULL)
#define P_INVALID (-1)

#define N_A 0
#define N_G 1
#define N_C 2
#define N_T 3
#define N_DUMMY 4

//Amino Acids Added
#define A_A (5)
#define A_R (6)
#define A_N (7)
#define A_D (8)
#define A_C (9)
#define A_E (10)
#define A_Q (11)
#define A_G (12)
#define A_H (13)
#define A_I (14)
#define A_L (15)
#define A_K (16)
#define A_M (17)
#define A_F (18)
#define A_P (19)
#define A_S (20)
#define A_T (21)
#define A_W (22)
#define A_Y (23)
#define A_V (24)
#define A_DUMMY (25)
#define A_INVALID (26)

#include "model.h"
#include "Utilities/LR_vector.h"
#include "Utilities/LR_matrix.h"
#include "Utilities/Logger.h"
#include "Utilities/parse_input/parse_input.h"

#include <string>
#include <memory>
#include <array>

using uint = uint32_t;
using llint = long long int;
using StressTensor = std::array<number, 6>;

#endif /* DEFS_H_ */
