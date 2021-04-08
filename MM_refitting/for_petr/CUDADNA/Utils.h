/*
 * Utils.h
 *
 *  Created on: 04/set/2010
 *      Author: lorenzo
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <cstdlib>
#include <cmath>
#include <cctype>

#include "defs.h"

class Utils {
public:
	Utils();
	virtual ~Utils();

	static int decode_base(char c);
	static char encode_base(int b);

	template<typename number> static number gaussian();
	template<typename number> static number sum(number *v, int N) {
		number res = (number) 0.;
		for(int i = 0; i < N; i++) res += v[i];
		return res;
	}

	template<typename number> static LR_vector<number> get_random_vector();
	template<typename number> static LR_vector<number> get_random_vector_in_sphere(number r);
	template<typename number> static void orthonormalize_matrix(LR_matrix<number> &M);
	template<typename number> static LR_matrix<number> get_random_rotation_matrix(number max_angle=2*M_PI);
	template<typename number> static LR_matrix<number> get_random_rotation_matrix_from_angle (number angle);
};

template<typename number>
inline LR_vector<number> Utils::get_random_vector() {
	number ransq = 1.;
	number ran1, ran2;

	while(ransq >= 1) {
		ran1 = 1. - 2. * drand48();
		ran2 = 1. - 2. * drand48();
		ransq = ran1*ran1 + ran2*ran2;
	}

	number ranh = 2. * sqrt(1. - ransq);
	return LR_vector<number>(ran1*ranh, ran2*ranh, 1. - 2. * ransq);
}

template<typename number>
inline LR_vector<number> Utils::get_random_vector_in_sphere(number r) {
	number r2 = SQR(r);
	LR_vector<number> res = LR_vector<number>(r, r, r);

	while(res.norm() > r2) {
		res = LR_vector<number>(2 * r * (drand48() - 0.5), 2 * r * (drand48() - 0.5), 2 * r * (drand48() - 0.5));
	}

	return res;
}

template<typename number>
void Utils::orthonormalize_matrix(LR_matrix<number> &m) {
    number v1_norm2 = m.v1 * m.v1;
    number v2_v1 = m.v2 * m.v1;

    m.v2 -= (v2_v1/v1_norm2) * m.v1;

    number v3_v1 = m.v3 * m.v1;
    number v3_v2 = m.v3 * m.v2;
    number v2_norm2 = m.v2 * m.v2;

    m.v3 -= (v3_v1/v1_norm2) * m.v1 + (v3_v2/v2_norm2) * m.v2;

    m.v1.normalize();
    m.v2.normalize();
    m.v3.normalize();
}

template<typename number>
LR_matrix<number> Utils::get_random_rotation_matrix_from_angle (number angle) {
	LR_vector<number> axis = Utils::get_random_vector<number>();

	number t = angle;
	number sintheta = sin(t);
	number costheta = cos(t);
	number olcos = 1. - costheta;

	number xyo = axis.x * axis.y * olcos;
	number xzo = axis.x * axis.z * olcos;
	number yzo = axis.y * axis.z * olcos;
	number xsin = axis.x * sintheta;
	number ysin = axis.y * sintheta;
	number zsin = axis.z * sintheta;

	LR_matrix<number> R(axis.x * axis.x * olcos + costheta, xyo - zsin, xzo + ysin,
				xyo + zsin, axis.y * axis.y * olcos + costheta, yzo - xsin,
				xzo - ysin, yzo + xsin, axis.z * axis.z * olcos + costheta);

//	LR_matrix<number> R(Utils::get_random_vector<number>(), Utils::get_random_vector<number>(), Utils::get_random_vector<number>());
//	Utils::orthonormalize_matrix<number>(R);
//
//	// rotations have det(R) == 1
//	if(R.determinant() < 0) {
//		LR_vector<number> tmp = R.v2;
//		R.v2 = R.v1;
//		R.v1 = tmp;
//	}

	return R;
}


template<typename number>
LR_matrix<number> Utils::get_random_rotation_matrix(number max_angle) {
	number t = max_angle * (drand48() - 0.5);
	return get_random_rotation_matrix_from_angle (t);
}

template<typename number>
inline number Utils::gaussian() {
	static unsigned int isNextG = 0;
	static number nextG;
	number toRet;
	number u, v, w;

	if(isNextG) {
		isNextG = 0;
		return nextG;
	}

	w = 2.;
	while(w >= 1.0) {
		u = 2. * drand48() - 1.0;
		v = 2. * drand48() - 1.0;
		w = u*u + v*v;
	}

	w = sqrt((-2. * log(w)) / w);
	toRet = u * w;
	nextG = v * w;
	isNextG = 1;

	return toRet;
}

#endif /* UTILS_H_ */
