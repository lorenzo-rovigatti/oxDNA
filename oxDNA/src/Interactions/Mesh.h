/*
 * Mesh; cubic mesh
 */

#ifndef MESH_H
#define MESH_H

/**
 * @brief Simple implementation of a cubic mesh.
 *
 * It is useful when dealing with interactions having computationally expensive
 * functional forms, such as Gaussian functions.
 */
template<typename number>
class Mesh {
public:
	Mesh() : N(0), A(NULL), B(NULL), C(NULL), D(NULL) {};

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
	number delta, inv_sqr_delta, xlow, xupp;
	number *A, *B, *C, *D;
};

#endif // MESH_H
