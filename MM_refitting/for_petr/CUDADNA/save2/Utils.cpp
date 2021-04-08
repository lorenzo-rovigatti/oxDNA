/*
 * Utils.cpp
 *
 *  Created on: 04/set/2010
 *      Author: lorenzo
 */

#include "Utils.h"

Utils::Utils() {

}

Utils::~Utils() {

}

int Utils::decode_base(char c) {
	c = toupper(c);
	switch(c) {
	case 'A': return N_A;
	case 'C': return N_C;
	case 'G': return N_G;
	case 'T': return N_T;
	case 'D': return N_DUMMY;
	default: return P_VIRTUAL;
	}
}

char Utils::encode_base(int b) {
	switch(b) {
	case N_A: return 'A';
	case N_C: return 'C';
	case N_G: return 'G';
	case N_T: return 'T';
	case N_DUMMY: return 'D';
	default: return 'X';
	}
}

template float Utils::gaussian<float>();
template double Utils::gaussian<double>();
