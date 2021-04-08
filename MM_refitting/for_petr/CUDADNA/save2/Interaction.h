/*
 * Interaction.h
 *
 *  Created on: Jun 7, 2011
 *      Author: rovigatti
 */

#ifndef INTERACTION_H_
#define INTERACTION_H_

#include <cassert>
#include <cstdio>
#include <cstdlib>

#include "defs.h"
extern "C" {
#include "parse_input/parse_input.h"
}

template<typename number>
class Mesh {
public:
     Mesh() : N(0) {};

	void init(int size) {
		N = size;
		delta = 0;
		A = new number[size + 1];
		B = new number[size + 1];
		C = new number[size + 1];
		D = new number[size + 1];
	}

	~Mesh() {
		if(N != 0) {
			delete[] A;
			delete[] B;
			delete[] C;
			delete[] D;
		}
	}

	int N;
	number delta, xlow, xupp;
	number *A, *B, *C, *D;
};

template <typename number>
class Interaction {
protected:
	bool _average;
	char _seq_filename[512];

	number _acosd(number x, void *a);
	number _acos(number x, void* a);
	void build_mesh(number (Interaction<number>::*f)(number, void*), number (Interaction<number>::*der) (number, void*), void * args, int npoints, number xlow, number xupp, Mesh<number> &m);

public:
	bool good;

	number F1_EPS[2][5][5];
	number F1_SHIFT[2][5][5];
	number F1_A[2];
	number F1_RC[2];
	number F1_R0[2];
	number F1_BLOW[2];
	number F1_BHIGH[2];
	number F1_RLOW[2];
	number F1_RHIGH[2];
	number F1_RCLOW[2];
	number F1_RCHIGH[2];

	number F2_K[2];
	number F2_RC[2];
	number F2_R0[2];
	number F2_BLOW[2];
	number F2_RLOW[2];
	number F2_RCLOW[2];
	number F2_BHIGH[2];
	number F2_RCHIGH[2];
	number F2_RHIGH[2];

	number F4_THETA_A[13];
	number F4_THETA_B[13];
	number F4_THETA_T0[13];
	number F4_THETA_TS[13];
	number F4_THETA_TC[13];

	number F5_PHI_A[4];
	number F5_PHI_B[4];
	number F5_PHI_XC[4];
	number F5_PHI_XS[4];

	number f1(number r, int type, int n3, int n5);
	number f1D(number r, int type, int n3, int n5);
	number f2(number r, int type);
	number f2D(number r, int type);
	number f4(number t, int type);
	number f4D(number t, int type);
	number f4DD(number t, int type);
	number f4Dsin(number t, int type);
	number fakef4(number t, void * par);
	number fakef4D(number t, void * par);
	number fakef4sp(number t, void * par);
	number fakef4Dsp(number t, void * par);
	number f5(number f, int type);
	number f5D(number f, int type);

	number lr_acos(number x);

	void get_settings(input_file &inp);
	void init(number T);

	number query_mesh(number x, Mesh<number> &m);
	number query_meshD(number x, Mesh<number> &m);

	int MESH_F4_POINTS[13];
	Mesh<number> mesh_f4[13];

	Interaction();
	virtual ~Interaction();
};

#endif /* INTERACTION_H_ */
