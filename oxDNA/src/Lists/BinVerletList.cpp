/*
 * BinVerletList.cpp
 *
 *  Created on: 07/nov/2013
 *      Author: lorenzo
 */

#include "BinVerletList.h"

using namespace std;

template<typename number>
BinVerletList<number>::BinVerletList(int &N, BaseBox<number> *box) : BaseList<number>(N, box), _updated(false) {

}

template<typename number>
BinVerletList<number>::~BinVerletList() {
	for(typename vector<Cells<number> *>::iterator it = _cells.begin(); it != _cells.end(); it++) delete *it;
}

template<typename number>
void BinVerletList<number>::get_settings(input_file &inp) {
	BaseList<number>::get_settings(inp);

	getInputNumber(&inp, "verlet_skin", &_skin, 1);
	_sqr_skin = SQR(_skin);

	getInputNumber(&inp, "bin_verlet_rcut[0]", _rcut, 1);
	getInputNumber(&inp, "bin_verlet_rcut[1]", _rcut + 1, 1);
	getInputNumber(&inp, "bin_verlet_rcut[2]", _rcut + 2, 1);

	if(this->_is_MC) {
		float delta_t = 0.f;
		getInputFloat(&inp, "delta_translation", &delta_t, 0);
		if(delta_t > 0.f && delta_t*sqrt(3) > _skin) throw oxDNAException("verlet_skin must be > delta_translation times sqrt(3) (the maximum displacement)");
	}
}

template<typename number>
void BinVerletList<number>::init(BaseParticle<number> **particles, number rcut) {
	this->_particles = particles;

	int N_part[2] = {0, 0};
	for(int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = particles[i];
		if(p->type != 0 && p->type != 1) throw oxDNAException("bin_verlet expects particles to be either of species 0 or 1, found %d", p->type);
		N_part[p->type]++;
	}

	if(N_part[0] == 0) OX_LOG(Logger::LOG_WARNING, "No particles of species 0 detected, why using bin_verlet then?");
	if(N_part[1] == 0) OX_LOG(Logger::LOG_WARNING, "No particles of species 1 detected, why using bin_verlet then?");

	for(int i = 0; i < 3; i++) {
		_rcut[i] += 2*_skin;
		_sqr_rcut[i] = SQR(_rcut[i]);
		Cells<number> *cn = new Cells<number>(this->_N, this->_box);
		cn->init(particles, _rcut[i]);
		_cells.push_back(cn);
	}

	_cells[0]->set_allowed_type(0);
	_cells[1]->set_unlike_type_only();
	_cells[2]->set_allowed_type(1);

	_lists.resize(this->_N, std::vector<BaseParticle<number> *>());
	_list_poss.resize(this->_N, LR_vector<number>(0, 0, 0));
	global_update(true);
}

template<typename number>
bool BinVerletList<number>::is_updated() {
	return (_updated);
}

template<typename number>
void BinVerletList<number>::single_update(BaseParticle<number> *p) {
	_cells[2*p->type]->single_update(p);
	_cells[1]->single_update(p);
	if(_list_poss[p->index].sqr_distance(p->pos) > _sqr_skin) _updated = false;
}

template<typename number>
void BinVerletList<number>::global_update(bool force_update) {
	if(!_cells[0]->is_updated() || !_cells[1]->is_updated() || !_cells[2]->is_updated() || force_update) {
		for(int i = 0; i < 3; i++) _cells[i]->global_update();
	}

	for(int i = 0; i < this->_N; i++) {
		BaseParticle<number> *p = this->_particles[i];
		_lists[p->index] = _cells[2*p->type]->get_neigh_list(p);
		vector<BaseParticle<number> *> l2 = _cells[1]->get_neigh_list(p);
		_lists[p->index].reserve(_lists[p->index].size() + l2.size());
		_lists[p->index].insert(_lists[p->index].end(), l2.begin(), l2.end());

		_list_poss[p->index] = p->pos;
	}
	_updated = true;
}

template<typename number>
std::vector<BaseParticle<number> *> BinVerletList<number>::get_neigh_list(BaseParticle<number> *p, bool all) {
	if(all) {
		vector<BaseParticle<number> *> l1 = _cells[2*p->type]->get_neigh_list(p, all);
		vector<BaseParticle<number> *> l2 = _cells[1]->get_neigh_list(p, all);

		l1.reserve(l1.size() + l2.size());
		l1.insert(l1.end(), l2.begin(), l2.end());
		return l1;
	}
	return _lists[p->index];
}

template class BinVerletList<float>;
template class BinVerletList<double>;
