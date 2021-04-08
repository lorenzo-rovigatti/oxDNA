/*
 * Weights.h
 *
 *  Created on: Feb 10, 2012
 *      Author: Flavio
 */

#ifndef WEIGHTS_H_
#define WEIGHTS_H_

#include "OrderParameters.h"

//using namespace std;

// Weight class
class Weights {
	protected:
		double * _w;
		int _dim;
		int _ndim;
		int * _sizes;
	public:
		Weights();
		~Weights();
		int get_dim (void) const { return _dim; }
		double get_weight_by_index (int);
		double get_weight(int *);
		double get_weight(int *, int *);
		double get_weight(OrderParameters *);
		double get_weight(OrderParameters *, int *);
		void init(const char *, OrderParameters *);
		void print();
};

#endif /* WEIGHTS_H_ */

