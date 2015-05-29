/*
 * VerletList.cpp
 *
 *  Created on: 07/nov/2013
 *      Author: lorenzo
 */

#include "VerletList.h"

template<typename number>
VerletList<number>::VerletList(int &N, BaseBox<number> *box) : BaseList<number>(N, box), _updated(false), _cells(N, box) {

}

template<typename number>
VerletList<number>::~VerletList() {

}

template<typename number>
void VerletList<number>::get_settings(input_file &inp) {
	BaseList<number>::get_settings(inp);
	_cells.get_settings(inp);

	getInputNumber(&inp, "verlet_skin", &_skin, 1);
	_sqr_skin = SQR(_skin);

	if(this->_is_MC) {
		float delta_t = 0.f;
		getInputFloat(&inp, "delta_translation", &delta_t, 0);
		if(delta_t > 0.f && delta_t * sqrt(3) > _skin) throw oxDNAException("verlet_skin must be > delta_translation times sqrt(3) (the maximum displacement)");
	}
}

template<typename number>
void VerletList<number>::init(BaseParticle<number> **particles, number rcut) {
	rcut += 2*_skin;
	BaseList<number>::init(particles, rcut);

	_sqr_rcut = SQR(rcut);

	_lists.resize(this->_N, std::vector<BaseParticle<number> *>());
	_list_poss.resize(this->_N, LR_vector<number>(0, 0, 0));

	_cells.init(particles, rcut);
	global_update();
}

template<typename number>
bool VerletList<number>::is_updated() {
	return (_updated);
}

template<typename number>
void VerletList<number>::single_update(BaseParticle<number> *p) {
	_cells.single_update(p);
	if(_list_poss[p->index].sqr_distance(p->pos) > _sqr_skin) _updated = false;
}

template<typename number>
void VerletList<number>::global_update(bool force_update) {
	if(!_cells.is_updated() || force_update) _cells.global_update();

	for(int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = this->_particles[i];
		_lists[p->index] = _cells.get_neigh_list(p);
		_list_poss[p->index] = p->pos;
	}
	_updated = true;
}

template<typename number>
std::vector<BaseParticle<number> *> VerletList<number>::get_neigh_list(BaseParticle<number> *p, bool all) {
	if(all) return _cells.get_neigh_list(p, true);
	return _lists[p->index];
}

template class VerletList<float>;
template class VerletList<double>;
