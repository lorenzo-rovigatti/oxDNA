/*
 * LR_matrix.h
 *
 *  Created on: 14 ott 2019
 *      Author: lorenzo
 */

#ifndef SRC_UTILITIES_LR_MATRIX_H_
#define SRC_UTILITIES_LR_MATRIX_H_

#include "LR_vector.h"

/**
 * @brief A 3x3 matrix with built-in basic matrix operations.
 *
 * This class, incapsulates a 3x3 matrix, internally stored as 3 {@link LR_vector}. It is not optimised
 * and hence should be either make it more efficient or substituted alltogether.
 */
class LR_matrix {
public:
	LR_vector v1, v2, v3;

	// Constructors
	LR_matrix(const LR_vector nv1, const LR_vector nv2, const LR_vector nv3) :
					v1(nv1),
					v2(nv2),
					v3(nv3) {
	}

	LR_matrix(const number n1, const number n2, const number n3, const number n4, const number n5, const number n6, const number n7, const number n8, const number n9) :
					v1(LR_vector(n1, n2, n3)),
					v2(LR_vector(n4, n5, n6)),
					v3(LR_vector(n7, n8, n9)) {
	}

	LR_matrix() {

	}

	LR_matrix(const LR_matrix &) = default;
	LR_matrix(LR_matrix &&) = default;
	LR_matrix &operator=(const LR_matrix &) = default;

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
		number fInv = ((number) 1.) / S;
		return LR_matrix(v1 * fInv, v2 * fInv, v3 * fInv);
	}

	inline void operator/=(number S) {
		v1 /= S;
		v2 /= S;
		v3 /= S;
	}

	inline LR_vector operator*(const LR_vector& v) const {
		return LR_vector(v1 * v, v2 * v, v3 * v);
	}

	inline LR_matrix operator*(number S) const {
		return LR_matrix(v1 * S, v2 * S, v3 * S);
	}

	inline LR_matrix operator*(const LR_matrix& m) const {
		LR_matrix tm = m.get_transpose();
		return LR_matrix(v1 * tm.v1, v1 * tm.v2, v1 * tm.v3, v2 * tm.v1, v2 * tm.v2, v2 * tm.v3, v3 * tm.v1, v3 * tm.v2, v3 * tm.v3);
	}

	inline LR_vector &operator[](int i) {
		if(i == 0) return v1;
		else if(i == 1) return v2;
		else return v3;
	}

	inline void transpone() {
		LR_vector nv1(v1.x, v2.x, v3.x);
		LR_vector nv2(v1.y, v2.y, v3.y);
		LR_vector nv3(v1.z, v2.z, v3.z);
		v1 = nv1;
		v2 = nv2;
		v3 = nv3;
	}

	inline LR_matrix get_transpose() const {
		return LR_matrix(v1.x, v2.x, v3.x, v1.y, v2.y, v3.y, v1.z, v2.z, v3.z);
	}

	inline number determinant() const {
		return v1.x * v2.y * v3.z + v1.y * v2.z * v3.x + v1.z * v2.x * v3.y - (v3.x * v2.y * v1.z + v3.y * v2.z * v1.x + v3.z * v2.x * v1.y);
	}

	inline void orthonormalize() {
		v1.normalize();
		v3.normalize();
		v1 -= v3 * (v1 * v3);
		v1.normalize();
		v2 = v3.cross(v1);
	}
	
	friend std::ostream &operator<<(std::ostream &os, LR_matrix const m) {
		return os << "[(" << m.v1.x << ", " << m.v1.y << ", " << m.v1.z << "),\n" << 
					 " (" << m.v2.x << ", " << m.v2.y << ", " << m.v2.z << "),\n" <<
					 " (" << m.v3.x << ", " << m.v3.y << ", " << m.v3.z << ")]";
	}
};

#endif /* SRC_UTILITIES_LR_MATRIX_H_ */
