/*
 * Mesh; cubic mesh
 */

#ifndef MESH_H
#define MESH_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvla"

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

	int N;
	number delta, inv_sqr_delta, xlow, xupp;
	std::vector<number> A, B, C, D;
};

#pragma GCC diagnostic pop

#endif // MESH_H
