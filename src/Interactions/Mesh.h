/*
 * Mesh; cubic mesh
 */

#ifndef MESH_H
#define MESH_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvla"

#include "../defs.h"

/**
 * @brief Simple implementation of a cubic mesh.
 *
 * It is useful when dealing with interactions having computationally expensive
 * functional forms, such as Gaussian functions.
 */

class Mesh {
public:
	Mesh() : N(0) {
		delta = inv_sqr_delta = xlow = xupp = -1.;
	};

	void init(int size) {
		N = size;
		delta = 0;
		A.resize(size + 1);
		B.resize(size + 1);
		C.resize(size + 1);
		D.resize(size + 1);
	}

	~Mesh() {

	}

	/**
	 * @brief Build a mesh by using a function and its derivative.
	 *
	 * @param f function
	 * @param der derivative of f
	 * @param pars pointer to a structure which contains function parameters
	 * @param npoints size of the mesh
	 * @param xlow the mesh is defined on a finite interval. This is the lower end
	 * @param xupp Upper end of the mesh interval
	 */
	void build(std::function<number(number, void*)> f, std::function<number(number, void*)> der, void *pars, int npoints, number xlow, number xupp);
	number query(number x);
	number query_derivative(number x);

	int N;
	number delta, inv_sqr_delta, xlow, xupp;
	std::vector<number> A, B, C, D;
};

#pragma GCC diagnostic pop

#endif // MESH_H
