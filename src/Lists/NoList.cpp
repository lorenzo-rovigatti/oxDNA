/*
 * NoList.cpp
 *
 *  Created on: 05/nov/2013
 *      Author: lorenzo
 */

#include "NoList.h"

template<typename number>
NoList<number>::NoList(int &N, number &box) : BaseList<number>(N, box) {

}

template<typename number>
NoList<number>::~NoList() {

}

template<typename number>
void NoList<number>::init(BaseParticle<number> **particles, number rcut) {
	BaseList<number>::init(particles, rcut);

	global_update();
}

template<typename number>
void NoList<number>::single_update(BaseParticle<number> *p) {

}

template<typename number>
void NoList<number>::global_update() {

}

template<typename number>
std::vector<BaseParticle<number> *> NoList<number>::get_neigh_list(BaseParticle<number> *p) {
	std::vector<BaseParticle<number> *> res;

	if(this->_is_MC) {
		for(int i = 0; i < this->_N; i++) if(i != p->index) res.push_back(this->_particles[i]);
	}
	else for(int i = p->index+1; i < this->_N; i++) res.push_back(this->_particles[i]);

	return res;
}

template class NoList<float>;
template class NoList<double>;
