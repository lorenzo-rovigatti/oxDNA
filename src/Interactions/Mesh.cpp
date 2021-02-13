/*
 * Mesh.cpp
 *
 *  Created on: Feb 13, 2021
 *      Author: lorenzo
 */

#include "Mesh.h"

#include <cfloat>
#include <cassert>

void Mesh::build(std::function<number(number, void*)> f, std::function<number(number, void*)> der, void *args, int npoints, number xlow, number xupp) {
	assert(xlow < xupp);
	int i;
	number x;

	this->init(npoints);

	number dx = (xupp - xlow) / (number) npoints;
	this->delta = dx;
	this->inv_sqr_delta = 1 / SQR(dx);
	this->xlow = xlow;
	this->xupp = xupp;

	number fx0, fx1, derx0, derx1;

	for(i = 0; i < npoints + 1; i++) {
		x = xlow + i * dx;

		fx0 = f(x, args);
		fx1 = f(x + dx, args);
		derx0 = der(x, args);
		derx1 = der(x + dx, args);

		this->A[i] = fx0;
		this->B[i] = derx0;
		this->D[i] = (2 * (fx0 - fx1) + (derx0 + derx1) * dx) / dx;
		this->C[i] = (fx1 - fx0 + (-derx0 - this->D[i]) * dx);
	}
}

number Mesh::query(number x) {
	if(x <= this->xlow) return this->A[0];
	if(x >= this->xupp) x = this->xupp - FLT_EPSILON;
	int i = (int) ((x - this->xlow) / this->delta);
	number dx = x - this->xlow - this->delta * i;
	return (this->A[i] + dx * (this->B[i] + dx * (this->C[i] + dx * this->D[i]) * this->inv_sqr_delta));
}

number Mesh::query_derivative(number x) {
	if(x < this->xlow) return this->B[0];
	if(x >= this->xupp) x = this->xupp - FLT_EPSILON;
	int i = (int) ((x - this->xlow) / this->delta);
	number dx = x - this->xlow - this->delta * i;
	return (this->B[i] + (2 * dx * this->C[i] + 3 * dx * dx * this->D[i]) * this->inv_sqr_delta);
}
