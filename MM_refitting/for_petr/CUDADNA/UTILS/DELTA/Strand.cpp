/*
 * Strand.cpp
 *
 *  Created on: 20/set/2011
 *      Author: lorenzo
 */

#include "Strand.h"

Strand::Strand(Particle<double> *first_particle, int length, int max_neighs) : cm_pos(LR_vector<double>(0, 0, 0)), n_neighs(0), neighs(NULL) {
	_particles = first_particle;
	_length = length;
	first = first_particle;
	last = _particles + length-1;
	_max_neighs = max_neighs;
	if(max_neighs > 0) neighs = new int[max_neighs];

	for(int i = 0; i < length; i++)	cm_pos += _particles[i].pos;

	cm_pos /= length;
}

Strand::~Strand() {
	if(neighs != NULL) delete[] neighs;
}

void Strand::add_neigh(int nn) {
	if(n_neighs < _max_neighs) {
		// check that nn is not inside the neighbour list already
		bool found = false;
		for(int i = 0; i < n_neighs; i++) if(neighs[i] == nn) found = true;
		if(!found) {
			neighs[n_neighs] = nn;
			n_neighs++;
		}
	}
	else {
		fprintf(stderr, "TROPPI VICINI!\n");
		exit(1);
	}
}

void Strand::translate(LR_vector<double> &disp) {
	for(int i = 0; i < _length; i++) {
		_particles[i].pos += disp;
		_particles[i].set_positions();
	}
	cm_pos += disp;
}

void Strand::rotate(LR_matrix<double> &R, LR_vector<double> &with_respect_to) {
	for(int i = 0; i < _length; i++) {
		_particles[i].pos = R * (_particles[i].pos - with_respect_to) + with_respect_to;
		_particles[i].orientation = R * _particles[i].orientation;
		_particles[i].orientationT = _particles[i].orientation.get_transpose();
		_particles[i].set_positions();
	}

	cm_pos = R*(cm_pos - with_respect_to) + with_respect_to;
}
