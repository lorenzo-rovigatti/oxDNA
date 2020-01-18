/*
 * BinVerletList.cpp
 *
 *  Created on: 07/nov/2013
 *      Author: lorenzo
 */

#include "BinVerletList.h"

using namespace std;

BinVerletList::BinVerletList(std::vector<BaseParticle *> &ps, BaseBox *box) :
				BaseList(ps, box),
				_updated(false),
				_is_AO(false) {

}

BinVerletList::~BinVerletList() {
	for(typename vector<Cells *>::iterator it = _cells.begin(); it != _cells.end(); it++)
		delete *it;
}

void BinVerletList::get_settings(input_file &inp) {
	BaseList::get_settings(inp);

	getInputNumber(&inp, "verlet_skin", &_skin, 1);
	_sqr_skin = SQR(_skin);

	getInputNumber(&inp, "bin_verlet_rcut[0]", _rcut, 1);
	getInputNumber(&inp, "bin_verlet_rcut[1]", _rcut + 1, 1);
	getInputNumber(&inp, "bin_verlet_rcut[2]", _rcut + 2, 1);

	getInputBool(&inp, "AO_mixture", &_is_AO, 0);

	if(this->_is_MC) {
		float delta_t = 0.f;
		getInputFloat(&inp, "delta_translation", &delta_t, 0);
		if(delta_t > 0.f && delta_t * sqrt(3) > _skin) throw oxDNAException("verlet_skin must be > delta_translation times sqrt(3) (the maximum displacement)");
	}
}

void BinVerletList::init(number rcut) {
	int N_part[2] = { 0, 0 };
	for(uint i = 0; i < _particles.size(); i++) {
		BaseParticle *p = _particles[i];
		if(p->type != 0 && p->type != 1) throw oxDNAException("bin_verlet expects particles to be either of species 0 or 1, found %d", p->type);
		N_part[p->type]++;
	}

	if(N_part[0] == 0) OX_LOG(Logger::LOG_WARNING, "No particles of species 0 detected, why using bin_verlet then?");
	if(N_part[1] == 0) OX_LOG(Logger::LOG_WARNING, "No particles of species 1 detected, why using bin_verlet then?");

	for(int i = 0; i < 3; i++) {
		_rcut[i] += 2 * _skin;
		_sqr_rcut[i] = SQR(_rcut[i]);
		Cells *cn = new Cells(_particles, _box);
		cn->init(_rcut[i]);
		_cells.push_back(cn);
	}

	_cells[0]->set_allowed_type(0);
	_cells[1]->set_unlike_type_only();
	_cells[2]->set_allowed_type(1);

	_lists.resize(_particles.size(), std::vector<BaseParticle *>());
	_list_poss.resize(_particles.size(), LR_vector(0, 0, 0));
	global_update(true);
}

bool BinVerletList::is_updated() {
	return (_updated);
}

void BinVerletList::single_update(BaseParticle *p) {
	_cells[2 * p->type]->single_update(p);
	_cells[1]->single_update(p);
	if(_list_poss[p->index].sqr_distance(p->pos) > _sqr_skin) _updated = false;
}

void BinVerletList::global_update(bool force_update) {
	if(!_cells[0]->is_updated() || !_cells[1]->is_updated() || !_cells[2]->is_updated() || force_update) {
		for(int i = 0; i < 3; i++)
			_cells[i]->global_update();
	}

	for(uint i = 0; i < _particles.size(); i++) {
		BaseParticle *p = this->_particles[i];
		if(p->type == 0 || !_is_AO) _lists[p->index] = _cells[2 * p->type]->get_neigh_list(p);
		else _lists[p->index].clear();
		vector<BaseParticle *> l2 = _cells[1]->get_neigh_list(p);
		_lists[p->index].reserve(_lists[p->index].size() + l2.size());
		_lists[p->index].insert(_lists[p->index].end(), l2.begin(), l2.end());

		_list_poss[p->index] = p->pos;
	}
	_updated = true;
}

std::vector<BaseParticle *> BinVerletList::get_neigh_list(BaseParticle *p) {
	return _lists[p->index];
}

std::vector<BaseParticle *> BinVerletList::get_complete_neigh_list(BaseParticle *p) {
	vector<BaseParticle *> l1 = _cells[2 * p->type]->get_complete_neigh_list(p);
	vector<BaseParticle *> l2 = _cells[1]->get_complete_neigh_list(p);

	l1.reserve(l1.size() + l2.size());
	l1.insert(l1.end(), l2.begin(), l2.end());
	return l1;
}
