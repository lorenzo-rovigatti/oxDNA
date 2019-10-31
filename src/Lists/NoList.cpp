/*
 * NoList.cpp
 *
 *  Created on: 05/nov/2013
 *      Author: lorenzo
 */

#include "NoList.h"


NoList::NoList(int &N, BaseBox *box) : BaseList(N, box) {

}


NoList::~NoList() {

}


void NoList::init(BaseParticle **particles, number rcut) {
	BaseList::init(particles, rcut);

	global_update(true);
}


void NoList::single_update(BaseParticle *p) {

}


void NoList::global_update(bool force_update) {

}


std::vector<BaseParticle *> NoList::_get_neigh_list(BaseParticle *p, bool all) {
	std::vector<BaseParticle *> res;

	int last = (all) ? this->_N : p->index;

	for(int i = 0; i < last; i++) {
		BaseParticle *q = this->_particles[i];
		if(p != q && !p->is_bonded(q)) res.push_back(this->_particles[i]);
	}

	return res;
}


std::vector<BaseParticle *> NoList::get_neigh_list(BaseParticle *p) {
	return _get_neigh_list(p, this->_is_MC);
}


std::vector<BaseParticle *> NoList::get_complete_neigh_list(BaseParticle *p) {
	return _get_neigh_list(p, true);
}
