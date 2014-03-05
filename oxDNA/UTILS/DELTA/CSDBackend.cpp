/*
 * CSDBackend.cpp
 *
 *  Created on: 26/set/2011
 *      Author: lorenzo
 */

#include "CSDBackend.h"
#include "../../IOManager.h"

CSDBackend::CSDBackend(IOManager *IO) : MC_CPUBackend<double>(IO) {

}

CSDBackend::~CSDBackend() {
	for(int i = 0; i < _N_cyls; i++) delete _cyls[i];

	delete[] _cyls;
	delete[] _clusters;
	delete[] _csd;
}

void CSDBackend::init(ifstream &conf_input) {
	MC_CPUBackend<double>::init(conf_input);

	if(_N % 24 != 0)_IO->die("The number of nucleotides must be a multiple of 24 (i.e. there must be only cylinders in the box)");

	_N_cyls = _N / 24;
	_cyls = new Strand*[_N_cyls];
	_clusters = new int[_N_cyls];
	_csd = new int[_N_cyls+1];

	for(int i = 0; i < _N_cyls; i++) {
		_cyls[i] = new Strand(_particles + i*24, 24, 4);
		_cyls[i]->cluster = i;
		_clusters[i] = 1;
		_csd[i] = 0;
	}
	_csd[_N_cyls] = 0;
}

void CSDBackend::_flip_neighs(int cyl) {
	Strand *c = _cyls[cyl];
	for(int i = 0; i < c->n_neighs; i++) {
		if(_cyls[c->neighs[i]]->cluster > c->cluster) {
			_clusters[_cyls[c->neighs[i]]->cluster]--;
			_cyls[c->neighs[i]]->cluster = c->cluster;
			_clusters[c->cluster]++;

			_flip_neighs(c->neighs[i]);
		}
	}
}

void CSDBackend::set_csd() {
	for(int i = 0; i < _N; i++) {
		Particle<double> *p = &_particles[i];
		p->prepare_list();
		int neigh = p->next_neighbour();
		while(neigh != P_VIRTUAL) {
			int c1 = i / 24;
			int c2 = neigh / 24;
			if(c2 > c1 && _particle_particle_interaction(p, &_particles[neigh]) < 0.) {
				_cyls[c1]->add_neigh(c2);
				_cyls[c2]->add_neigh(c1);
			}
			neigh = p->next_neighbour();
		}
	}

	for(int i = 0; i < _N_cyls; i++) _flip_neighs(i);
	for(int i = 0; i < _N_cyls; i++) _csd[_clusters[i]]++;
}
