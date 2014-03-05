/**
 * @file    Weights.h
 * @date    Feb 10, 2012
 * @author  flavio
 * 
 */

#ifndef WEIGHTS_H_
#define WEIGHTS_H_

#include "OrderParameters.h"

/// Weight class
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
	void init(const char *, OrderParameters *, bool safe, double default_weight);
	void print();
};

#endif /* WEIGHTS_H_ */

