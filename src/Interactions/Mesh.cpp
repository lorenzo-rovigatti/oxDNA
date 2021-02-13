/*
 * Mesh.cpp
 *
 *  Created on: Feb 13, 2021
 *      Author: lorenzo
 */

#include "Mesh.h"

#include <cfloat>
#include <cassert>
#include <functional>

void Mesh::build(std::function<number(number, void*)> f, std::function<number(number, void*)> der, void *args, int npoints, number xlow, number xupp) {
	assert(xlow < xupp);
	int i;
	number x;

	this->init(npoints);

	number dx = (xupp - xlow) / (number) npoints;
	this->_delta = dx;
	this->_inv_sqr_delta = 1 / SQR(dx);
	this->_xlow = xlow;
	this->_xupp = xupp;

	number fx0, fx1, derx0, derx1;

	for(i = 0; i < npoints + 1; i++) {
		x = xlow + i * dx;

		fx0 = f(x, args);
		fx1 = f(x + dx, args);
		derx0 = der(x, args);
		derx1 = der(x + dx, args);

		this->_A[i] = fx0;
		this->_B[i] = derx0;
		this->_D[i] = (2 * (fx0 - fx1) + (derx0 + derx1) * dx) / dx;
		this->_C[i] = (fx1 - fx0 + (-derx0 - this->_D[i]) * dx);
	}
}

number Mesh::query(number x) {
	if(x <= this->_xlow) return this->_A[0];
	if(x >= this->_xupp) x = this->_xupp - FLT_EPSILON;
	int i = (int) ((x - this->_xlow) / this->_delta);
	number dx = x - this->_xlow - this->_delta * i;
	return (this->_A[i] + dx * (this->_B[i] + dx * (this->_C[i] + dx * this->_D[i]) * this->_inv_sqr_delta));
}

number Mesh::query_derivative(number x) {
	if(x < this->_xlow) return this->_B[0];
	if(x >= this->_xupp) x = this->_xupp - FLT_EPSILON;
	int i = (int) ((x - this->_xlow) / this->_delta);
	number dx = x - this->_xlow - this->_delta * i;
	return (this->_B[i] + (2 * dx * this->_C[i] + 3 * dx * dx * this->_D[i]) * this->_inv_sqr_delta);
}
