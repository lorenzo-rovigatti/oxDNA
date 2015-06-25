/*
 * NoList.cpp
 *
 *  Created on: 05/nov/2013
 *      Author: lorenzo
 */

#include "NoList.h"

template<typename number>
NoList<number>::NoList(int &N, BaseBox<number> *box) : BaseList<number>(N, box) {

}

template<typename number>
NoList<number>::~NoList() {

}

template<typename number>
void NoList<number>::init(BaseParticle<number> **particles, number rcut) {
	BaseList<number>::init(particles, rcut);

	global_update(true);
}

template<typename number>
void NoList<number>::single_update(BaseParticle<number> *p) {

}

template<typename number>
void NoList<number>::global_update(bool force_update) {

}

template<typename number>
std::vector<BaseParticle<number> *> NoList<number>::get_neigh_list(BaseParticle<number> *p, bool all) {
	std::vector<BaseParticle<number> *> res;

	int last = p->index;
	if(this->_is_MC && !all) last = this->_N;

	for(int i = 0; i < last; i++) {
		BaseParticle<number> *q = this->_particles[i];
		if(p != q && !p->is_bonded(q)) res.push_back(this->_particles[i]);
	}

	return res;
}

template class NoList<float>;
template class NoList<double>;
