/**
 * @file    defs.h
 * @date    13/ott/2009
 * @author  lorenzo
 *
 *
 */

#ifndef DEFS_H_
#define DEFS_H_

#define VERSION_MAJOR 2
#define VERSION_MINOR 3
#define VERSION_STAGE 2 // 0 alpha, 1 beta, 2 stable

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

#include <cmath>
#include <strings.h>
#include <iostream>

#include "model.h"
#include "Utilities/Logger.h"
#include "Utilities/parse_input/parse_input.h"

typedef long long int llint;
template<typename number> class LR_vector;

/**
 * @brief A 3x3 matrix with built-in basic matrix operations.
 *
 * This class, incapsulates a 3x3 matrix, internally stored as 3 {@link LR_vector}. It is not optimised
 * and hence should be either make it more efficient or substituted alltogether.
 */
template<typename number> class LR_matrix {
public:
	LR_vector<number> v1, v2, v3;

	// Constructors
	LR_matrix(const LR_vector<number> nv1, const LR_vector<number> nv2, const LR_vector<number> nv3) :
		v1(nv1), v2(nv2), v3(nv3) {
	}

	LR_matrix(const number n1, const number n2, const number n3, const number n4, const number n5,
			const number n6, const number n7, const number n8, const number n9) :
		v1(LR_vector<number>(n1, n2, n3)), v2(LR_vector<number>(n4, n5, n6)), v3(LR_vector<number>(n7, n8, n9)) {
	}

	LR_matrix() {

	}

	// Operator Overloads
	inline LR_matrix operator+(const LR_matrix& m) const {
		return LR_matrix(v1 + m.v1, v2 + m.v2, v3 + m.v3);
	}

	inline LR_matrix operator-(const LR_matrix& m) const {
		return LR_matrix(v1 - m.v1, v2 - m.v2, v3 - m.v3);
	}

	inline LR_matrix operator-() const {
		return LR_matrix(-v1, -v2, -v3);
	}

	inline LR_matrix operator+() const {
		return LR_matrix(v1, v2, v3);
	}

	inline LR_matrix operator/(number S) const {
		number fInv = ((number)1.) / S;
		return LR_matrix(v1 * fInv, v2 * fInv, v3 * fInv);
	}

	inline void operator/=(number S) {
			v1 /= S;
			v2 /= S;
			v3 /= S;
		}

	inline LR_vector<number> operator*(const LR_vector<number>& v) const {
		return LR_vector<number>(v1 * v, v2 * v, v3 * v);
	}

	inline LR_matrix operator*(number S) const {
		return LR_matrix(v1 * S, v2 * S, v3 * S);
	}

	inline LR_matrix operator*(const LR_matrix& m) const {
		LR_matrix tm = m.get_transpose();
		return LR_matrix<number>(
				v1 * tm.v1, v1 * tm.v2, v1 * tm.v3,
				v2 * tm.v1, v2 * tm.v2, v2 * tm.v3,
				v3 * tm.v1, v3 * tm.v2, v3 * tm.v3);
	}
	
	inline LR_vector<number> &operator[](int i) {
		if (i==0) return v1;
		else if (i==1) return v2;
		else return v3;
	}

	inline void transpone() {
		LR_vector<number> nv1 = LR_vector<number>(v1.x, v2.x, v3.x);
		LR_vector<number> nv2 = LR_vector<number>(v1.y, v2.y, v3.y);
		LR_vector<number> nv3 = LR_vector<number>(v1.z, v2.z, v3.z);
		v1 = nv1;
		v2 = nv2;
		v3 = nv3;
	}

	inline LR_matrix get_transpose() const {
		return LR_matrix(v1.x, v2.x, v3.x, v1.y, v2.y, v3.y, v1.z, v2.z, v3.z);
	}

	inline number determinant() const {
		return v1.x*v2.y*v3.z + v1.y*v2.z*v3.x + v1.z*v2.x*v3.y - (v3.x*v2.y*v1.z + v3.y*v2.z*v1.x + v3.z*v2.x*v1.y);
	}

	inline void orthonormalize() {
		v1.normalize ();
		v3.normalize ();
		v1 -= v3 * (v1 * v3);
		v1.normalize ();
		v2 = v3.cross (v1);
	}
};

template <typename number>
LR_vector<number> operator*(const number S, const LR_vector<number> &v);

/**
 * @brief A three-dimensional vector with built-in basic vectorial operations.
 */
template<typename number> class LR_vector {
public:
	// Data
	number x, y, z;

	// Constructors
	LR_vector(number InX, number InY, number InZ) :
		x(InX), y(InY), z(InZ) {
	}

	LR_vector() :
		x(0), y(0), z(0) {
	}

	LR_vector (number * arg) : 
		x (arg[0]), y (arg[1]), z (arg[2]) {
	
	}

	// Operator Overloads
	inline LR_vector operator+(const LR_vector& V2) const {
		return LR_vector(x + V2.x, y + V2.y, z + V2.z);
	}

	inline LR_vector operator-(const LR_vector& V2) const {
		return LR_vector(x - V2.x, y - V2.y, z - V2.z);
	}

	inline LR_vector operator-() const {
		return LR_vector(-x, -y, -z);
	}

	inline LR_vector operator+() const {
		return LR_vector(x, y, z);
	}

	inline LR_vector operator/(number S) const {
		number fInv = ((number)1.) / S;
		return LR_vector(x * fInv, y * fInv, z * fInv);
	}

	inline number operator*(const LR_vector& V2) const {
		return x * V2.x + y * V2.y + z * V2.z;
	}

	inline LR_vector operator*(number S) const {
		return LR_vector(x * S, y * S, z * S);
	}

	inline void operator+=(const LR_vector& V2) {
		x += V2.x;
		y += V2.y;
		z += V2.z;
	}
	inline void operator-=(const LR_vector& V2) {
		x -= V2.x;
		y -= V2.y;
		z -= V2.z;
	}

	inline void operator/=(number S) {
		x /= S;
		y /= S;
		z /= S;
	}

	inline void operator*=(number S) {
		x *= S;
		y *= S;
		z *= S;
	}

	inline number &operator[](int i) {
		if (i == 0)
			return x;
		else if (i == 1)
			return y;
		else
			return z;
	}

	inline LR_vector cross(const LR_vector &V2) const {
		return LR_vector(y * V2.z - z * V2.y, z * V2.x - x * V2.z, x * V2.y - y	* V2.x);
	}

	number norm() const {
		return x * x + y * y + z * z;
	}

	number module() const {
		return sqrt(norm());
	}

	number sqr_distance(const LR_vector &V1) const {
		return (*this - V1).norm();
	}

	number distance(const LR_vector &V1) const {
		return sqrt(sqr_distance(V1));
	}

	inline void normalize() {
		number fMag = (x * x + y * y + z * z);
		if (fMag == 0) return;

		number fMult = ((number)1.) / sqrtf(fMag);
		x *= fMult;
		y *= fMult;
		z *= fMult;
		return;
	}

	friend std::ostream &operator<<(std::ostream &stream,const  LR_vector &vector) {
		stream << "(" << vector.x << ", " << vector.y << ", " << vector.z << ")";
		return stream;
	}
};

template <typename number>
LR_vector<number> operator*(const number S, const LR_vector<number> &v) {
	return v*S;
}

#endif /* DEFS_H_ */
