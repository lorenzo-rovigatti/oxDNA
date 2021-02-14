/*
 * Mesh; cubic mesh
 */

#ifndef MESH_H
#define MESH_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvla"

#include "../defs.h"

#include <functional>
#include <cfloat>

/**
 * @brief Simple implementation of a cubic mesh.
 *
 * It is useful when dealing with interactions having computationally expensive
 * functional forms, such as Gaussian functions.
 */
class Mesh {
public:
	Mesh() : _N(0) {
		_delta = _inv_sqr_delta = _xlow = _xupp = -1.;
	};

	void init(int size) {
		_N = size;
		_delta = 0;
		_A.resize(size + 1);
		_B.resize(size + 1);
		_C.resize(size + 1);
		_D.resize(size + 1);
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

	inline number query(number x) {
		if(x <= _xlow) return _A[0];
		if(x >= _xupp) x = _xupp - FLT_EPSILON;
		int i = (int) ((x - _xlow) / _delta);
		number dx = x - _xlow - _delta * i;
		return (_A[i] + dx * (_B[i] + dx * (_C[i] + dx * _D[i]) * _inv_sqr_delta));
	}

	inline number query_derivative(number x) {
		if(x < _xlow) return _B[0];
		if(x >= _xupp) x = _xupp - FLT_EPSILON;
		int i = (int) ((x - _xlow) / _delta);
		number dx = x - _xlow - _delta * i;
		return (_B[i] + (2 * dx * _C[i] + 3 * dx * dx * _D[i]) * _inv_sqr_delta);
	}

	number x_low() {
		return _xlow;
	}

private:
	int _N;
	number _delta, _inv_sqr_delta, _xlow, _xupp;
	std::vector<number> _A, _B, _C, _D;
};

#pragma GCC diagnostic pop

#endif // MESH_H
