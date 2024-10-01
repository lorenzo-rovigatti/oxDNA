/*
 * Mesh.cpp
 *
 *  Created on: Feb 13, 2021
 *      Author: lorenzo
 */

#include "Mesh.h"

#include <cassert>

void Mesh::build(std::function<number(number, void*)> f, std::function<number(number, void*)> der, void *args, int npoints, number xlow, number xupp) {
	assert(xlow < xupp);
	int i;
	number x;

	init(npoints);

	number dx = (xupp - xlow) / (number) npoints;
	_delta = dx;
	_inv_sqr_delta = 1 / SQR(dx);
	_xlow = xlow;
	_xupp = xupp;

	number fx0, fx1, derx0, derx1;

	for(i = 0; i < npoints + 1; i++) {
		x = xlow + i * dx;

		fx0 = f(x, args);
		fx1 = f(x + dx, args);
		derx0 = der(x, args);
		derx1 = der(x + dx, args);

		_A[i] = fx0;
		_B[i] = derx0;
		_D[i] = (2 * (fx0 - fx1) + (derx0 + derx1) * dx) / dx;
		_C[i] = (fx1 - fx0 + (-derx0 - _D[i]) * dx);
	}
}
