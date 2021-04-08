/*
 * Strand.h
 *
 *  Created on: 20/set/2011
 *      Author: lorenzo
 */

#ifndef STRAND_H_
#define STRAND_H_

#include "../../Particle.h"

class Strand {
protected:
	Particle<double> *_particles;
	int _length;
	int _max_neighs;

public:
	LR_vector<double> cm_pos;
	Particle<double> *first;
	Particle<double> *last;

	int cluster;
	int n_neighs;
	int *neighs;

	Strand(Particle<double> *first_particle, int length, int max_neighs=0);
	virtual ~Strand();

	virtual void add_neigh(int nn);
	virtual void translate(LR_vector<double> &disp);
	virtual void rotate(LR_matrix<double> &R, LR_vector<double> &with_respect_to);
	virtual void rotate(LR_matrix<double> &R) { rotate(R, cm_pos); }
};

#endif /* STRAND_H_ */
