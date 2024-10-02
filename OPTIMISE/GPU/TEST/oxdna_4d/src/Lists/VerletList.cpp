/*
 * VerletList.cpp
 *
 *  Created on: 07/nov/2013
 *      Author: lorenzo
 */

#include "VerletList.h"

VerletList::VerletList(std::vector<BaseParticle *> &ps, BaseBox *box) :
				BaseList(ps, box),
				_updated(false),
				_cells(ps, box) {

}

VerletList::~VerletList() {

}

void VerletList::get_settings(input_file &inp) {
	BaseList::get_settings(inp);
	_cells.get_settings(inp);

	getInputNumber(&inp, "verlet_skin", &_skin, 1);
	_sqr_skin = SQR(_skin);

	if(this->_is_MC) {
		float delta_t = 0.f;
		getInputFloat(&inp, "delta_translation", &delta_t, 0);
		if(delta_t > 0.f && delta_t * sqrt(3) > _skin) throw oxDNAException("verlet_skin must be > delta_translation times sqrt(3) (the maximum displacement)");
	}
}

void VerletList::init(number rcut) {
	rcut += 2 * _skin;
	BaseList::init(rcut);

	_sqr_rcut = SQR(rcut);

	_lists.resize(_particles.size(), std::vector<BaseParticle *>());
	_list_poss.resize(_particles.size(), LR_vector(0, 0, 0));

	_cells.init(rcut);
	global_update();
}

bool VerletList::is_updated() {
	return (_updated);
}

void VerletList::single_update(BaseParticle *p) {
	_cells.single_update(p);
	if(_list_poss[p->index].sqr_distance(p->pos) > _sqr_skin) _updated = false;
}

void VerletList::global_update(bool force_update) {
	if(!_cells.is_updated() || force_update) _cells.global_update();

	for(uint i = 0; i < _particles.size(); i++) {
		BaseParticle *p = this->_particles[i];
		_lists[p->index] = _cells.get_neigh_list(p);
		_list_poss[p->index] = p->pos;
	}
	_updated = true;
}

std::vector<BaseParticle *> VerletList::get_neigh_list(BaseParticle *p) {
	return _lists[p->index];
}

std::vector<BaseParticle *> VerletList::get_complete_neigh_list(BaseParticle *p) {
	return _cells.get_complete_neigh_list(p);
}

void VerletList::change_box() {
	LR_vector new_box_sides = this->_box->box_sides();
	number fx = new_box_sides.x / this->_box_sides.x;
	number fy = new_box_sides.y / this->_box_sides.y;
	number fz = new_box_sides.z / this->_box_sides.z;

	for(uint i = 0; i < _particles.size(); i++) {
		BaseParticle *p = this->_particles[i];
		_list_poss[p->index].x *= fx;
		_list_poss[p->index].z *= fy;
		_list_poss[p->index].y *= fz;

		if(_list_poss[p->index].sqr_distance(p->pos) > _sqr_skin) _updated = false;
	}

	_cells.change_box();
	BaseList::change_box();
}
